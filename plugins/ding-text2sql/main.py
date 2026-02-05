import hashlib
import json
from typing import Any

from fastmcp import Context, FastMCP
from sqlalchemy import create_engine, text, inspect
from sqlalchemy.engine import Engine
from fastmcp.dependencies import CurrentContext, Depends


mcp = FastMCP("Dynamic Text2SQL Service")


class ConnectionManager:
    """
    负责管理数据库连接引擎的缓存。
    """

    def __init__(self):
        self._engines: dict[str, Engine] = {}

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


def get_config(ctx: Context = CurrentContext()) -> dict[str, Any]:
    meta = ctx.request_context.meta
    assert meta is not None
    plugin_config = meta.model_dump()
    return plugin_config


class DisplayConfig:
    """
    定义数据展示配置
    """

    @staticmethod
    def table(
        columns: list[str] | None = None, title: str | None = None
    ) -> dict[str, Any]:
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


def build_table_display(
    rows: list[dict[str, Any]], display_config: dict[str, Any] | None = None
) -> dict[str, Any] | None:
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


@mcp.tool()
def list_tables(
    plugin_config=Depends(get_config),
):
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


from sqlalchemy import inspect
from sqlalchemy.types import Integer, Numeric, String, Text, Date, DateTime, Boolean


@mcp.tool()
def get_database_schema(
    table_names: list[str] | None = None,
    plugin_config=Depends(get_config),
) -> str:
    """
    获取极致压缩的数据库 Schema，专为 LLM 优化。
    格式: TableName(col:type*PK, col:type>FK_Target, ...)
    """
    try:
        engine = manager.get_engine(plugin_config)
        inspector = inspect(engine)
        all_tables = inspector.get_table_names()

        # 如果未指定表名，只返回简单的表列表，不加载列
        if not table_names:
            return f"Available Tables: {', '.join(all_tables)}"

        schema_lines = []

        # 简单的类型映射字典，将复杂的 SQL 类型映射为 Python/TS 基础类型
        # LLM 只需要知道是 'int' 还是 'str' 就能写 SQL，不需要知道是 VARCHAR(255)
        type_map = {
            Integer: "int",
            Numeric: "float",
            String: "str",
            Text: "str",
            Date: "date",
            DateTime: "datetime",
            Boolean: "bool",
        }

        for table in table_names:
            if table not in all_tables:
                continue

            # 1. 获取外键映射
            try:
                fks = inspector.get_foreign_keys(table)
                fk_map = {}
                for fk in fks:
                    # 简化外键描述：只保留目标表名.列名
                    if (
                        fk.get("constrained_columns")
                        and fk.get("referred_table")
                        and fk.get("referred_columns")
                    ):
                        src_col = fk["constrained_columns"][0]
                        target = f"{fk['referred_table']}.{fk['referred_columns'][0]}"
                        fk_map[src_col] = target
            except:
                fk_map = {}

            # 2. 处理列信息
            columns = inspector.get_columns(table)
            col_strs = []

            for col in columns:
                col_name = col["name"]

                # 智能简化类型名称
                col_type_cls = type(col["type"])
                # 默认取类名的小写 (例如 INTEGER -> integer)
                simple_type = col_type_cls.__name__.lower()
                # 如果在映射表中，则使用更短的别名
                for base_type, alias in type_map.items():
                    if issubclass(col_type_cls, base_type):
                        simple_type = alias
                        break

                # 移除 'varchar' 等类型中可能包含的 '()' 防止 WAF 误判
                simple_type = simple_type.split("(")[0]

                # 组装属性：name:type
                part = f"{col_name}:{simple_type}"

                # 添加主键标记 (*PK)
                if col.get("primary_key"):
                    part += "*PK"

                # 添加外键标记 (>Target)
                if col_name in fk_map:
                    part += f">{fk_map[col_name]}"

                col_strs.append(part)

            # 3. 组装单行 Schema
            # 格式: table_name(col1:type, col2:type...)
            schema_lines.append(f"{table}({', '.join(col_strs)})")

        return "\n".join(schema_lines)

    except Exception:
        return "Error: Unable to retrieve schema."


@mcp.tool()
def execute_sql(
    ctx: Context,
    query: str,
    display: dict[str, Any] | None = None,
    plugin_config=Depends(get_config),
) -> str | dict[str, Any]:
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
    meta = ctx.request_context.meta

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
