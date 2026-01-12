import hashlib
import json
from typing import Dict, List, Optional

from fastmcp import FastMCP
from sqlalchemy import create_engine, text, inspect
from sqlalchemy.engine import Engine

mcp = FastMCP("Dynamic Text2SQL Service")


class ConnectionManager:
    """
    负责管理数据库连接引擎的缓存。
    """

    def __init__(self):
        self._engines: Dict[str, Engine] = {}

    def _get_cache_key(self, config: dict[str, str]) -> str:
        """
        根据字典内容生成唯一指纹。
        要求 config 字典中至少包含 'db_uri'。
        """
        config_str = json.dumps(config, sort_keys=True)
        return hashlib.md5(config_str.encode()).hexdigest()

    def get_engine(self, config: dict[str, str]) -> Engine:
        conn_str = config.get("db_uri")

        if not conn_str:
            raise ValueError("Config missing required key:  'db_uri'")

        if conn_str.strip().lower().startswith("sqlite"):
            raise ValueError(
                "SQLite databases are currently disabled based on security policy."
            )

        key = self._get_cache_key(config)
        if key in self._engines:
            return self._engines[key]

        engine = create_engine(conn_str, pool_recycle=3600)
        self._engines[key] = engine
        return engine


# 初始化全局连接管理器
manager = ConnectionManager()


# --- 展示配置类 ---


class DisplayConfig:
    """
    定义数据展示配置
    """

    @staticmethod
    def table(columns: Optional[List[str]] = None, title: Optional[str] = None) -> dict:
        """
        创建表格展示配置

        Args:
            columns: 要展示的列名列表，None 表示展示所有列
            title: 表格标题（可选）
        """
        config = {"type": "table"}
        if columns:
            config["columns"] = columns
        if title:
            config["title"] = title
        return config


from typing import List, Optional, Dict, Any


def build_table_display(
    rows: List[dict], display_config: Optional[dict] = None
) -> dict:
    """
    根据配置构建表格展示数据 (适配 accessorKey 前端模式)

    Args:
        rows: 查询结果行列表 (List[dict])
        display_config: 展示配置

    Returns:
        表格展示格式的字典，包含 columns 配置和 rows 数据
    """
    if not rows:
        return None

    # 1. 确定要展示的列名 (target_columns)
    all_columns = list(rows[0].keys())
    target_columns = all_columns

    if display_config and display_config.get("columns"):
        # 过滤出存在于数据中的列，并保持配置顺序
        target_columns = [
            col for col in display_config["columns"] if col in all_columns
        ]
        # 如果配置的列都不存在，回退到展示所有列
        if not target_columns:
            target_columns = all_columns

    # 2. 构建列定义 (Column Definitions)
    # 对应你之前给出的结构: [{"accessorKey": "breed_cn_name"}, ...]

    # 3. 构建行数据 (Data Rows)
    # 关键修改：保持为对象列表 List[dict]，而非 List[list]
    table_rows = []
    for row in rows:
        # 只保留需要的键值对
        filtered_row = {col: row.get(col) for col in target_columns}
        table_rows.append(filtered_row)

    result = {
        "type": "table",
        "columns": target_columns,  # 前端拿去渲染表头
        "rows": table_rows,  # 前端拿去填充数据
    }

    # 添加可选的标题
    if display_config and display_config.get("title"):
        result["title"] = display_config["title"]

    return result


# --- MCP 工具定义 ---


@mcp.tool(exclude_args=["plugin_config"])
def list_tables(plugin_config: dict[str, str] = {}) -> str:
    """
    在使用 `get_table_schema` 之前，先使用此工具查看有哪些表。
    """
    try:
        engine = manager.get_engine(plugin_config)
        inspector = inspect(engine)
        table_names = inspector.get_table_names()

        if not table_names:
            return "Database is empty (no tables found)."

        return json.dumps(table_names, ensure_ascii=False)

    except Exception as e:
        return f"Error listing tables: {str(e)}"


@mcp.tool(exclude_args=["plugin_config"])
def get_database_schema(
    table_names: list[str] | None = None,
    plugin_config: dict[str, str] = {},
) -> str:
    """
    获取数据库 Schema。
    """
    try:
        engine = manager.get_engine(plugin_config)
        inspector = inspect(engine)
        all_tables = inspector.get_table_names()

        if not table_names:
            return f"Found tables: {', '.join(all_tables)}\n\nPlease call this tool again with specific `table_names` to get column details."

        schema_info = []
        for table in table_names:
            if table not in all_tables:
                schema_info.append(f"Table '{table}' not found.")
                continue

            columns = inspector.get_columns(table)
            col_desc = []
            for col in columns:
                col_str = f"- {col['name']} ({str(col['type'])})"
                if col.get("primary_key"):
                    col_str += " [PK]"
                if col.get("foreign_keys"):
                    fk_texts = [
                        f"FK -> {fk.column.table.name}.{fk.column.name}"
                        for fk in col["foreign_keys"]
                    ]
                    col_str += f" [{', '.join(fk_texts)}]"
                col_desc.append(col_str)

            schema_info.append(f"Table: {table}\n" + "\n".join(col_desc))

        return "\n\n".join(schema_info)

    except Exception as e:
        return f"Error getting schema: {str(e)}"


@mcp.tool(exclude_args=["plugin_config"])
def execute_sql(
    query: str, display: Optional[dict] = None, plugin_config: dict[str, str] = {}
) -> str | dict:
    """
    执行 SQL 查询并返回结果。

    Args:
        query:  要执行的 SQL Select 语句。
        display: 可视化展示配置，格式如下：
            - {"type": "table"}:  以表格形式展示所有列
            - {"type": "table", "columns": ["col1", "col2"]}: 以表格形式展示指定列
            - {"type":  "table", "title": "标题"}: 带标题的表格
            - None: 不进行可视化展示，只返回原始数据

    Returns:
        返回格式为 {"model_text": text, "display": display_list}
        其中 display_list 是展示配置列表
    """
    normalized_query = query.strip().lower()

    # 基础安全检查
    forbidden = [
        "drop ",
        "delete ",
        "update ",
        "insert ",
        "alter ",
        "grant ",
        "truncate ",
    ]
    if any(cmd in normalized_query for cmd in forbidden):
        return "Error: Safety restriction - Only SELECT queries are allowed."

    try:
        engine = manager.get_engine(plugin_config)
        with engine.connect() as conn:
            result = conn.execute(text(query))
            keys = list(result.keys())
            rows = result.fetchmany(50)  # 限制返回行数

            if not rows:
                return {"model_text": "Query returned no results.", "display": []}

            # 转换为字典列表
            output = [dict(zip(keys, row)) for row in rows]
            model_text = json.dumps(output, default=str, ensure_ascii=False)

            # 构建展示数据
            display_list = []

            if display:
                display_type = display.get("type", "table")

                if display_type == "table":
                    table_display = build_table_display(output, display)
                    if table_display:
                        display_list.append(table_display)
                # 未来可以在这里扩展其他类型，如 chart, graph 等
                # elif display_type == "chart":
                #     chart_display = build_chart_display(output, display)
                #     if chart_display:
                #         display_list.append(chart_display)

            return {"model_text": model_text, "display": display_list}

    except Exception as e:
        return f"SQL Execution Error: {str(e)}"


if __name__ == "__main__":
    mcp.run()
