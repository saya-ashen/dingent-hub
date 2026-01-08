import hashlib
import json
from typing import Dict, List

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
        # 为了保证指纹一致性，我们需要对字典 key 进行排序
        config_str = json.dumps(config, sort_keys=True)
        return hashlib.md5(config_str.encode()).hexdigest()

    def get_engine(self, config: dict[str, str]) -> Engine:
        conn_str = config.get("db_uri")

        if not conn_str:
            raise ValueError("Config missing required key: 'db_uri'")
        # ------------------- 限制逻辑 -------------------
        # 检查是否为 SQLite 连接字符串
        # 通常格式为 sqlite:///:memory: 或 sqlite:///relative/path.db 或 sqlite:////absolute/path.db
        if conn_str.strip().lower().startswith("sqlite"):
            # TODO: 未来如果需要开放固定路径，可以在此处解析 path 并进行白名单校验
            # 例如：
            # if "/allowed/path/db" not in conn_str:
            #     raise ValueError("SQLite path not allowed")

            raise ValueError(
                "SQLite databases are currently disabled based on security policy."
            )
        # ----------------------------------------------------

        # 2. 检查缓存
        key = self._get_cache_key(config)
        if key in self._engines:
            return self._engines[key]

        engine = create_engine(conn_str, pool_recycle=3600)
        self._engines[key] = engine
        return engine


# 初始化全局连接管理器
manager = ConnectionManager()

# --- 3. MCP 工具定义 ---


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

        # 返回 JSON 格式的列表，LLM 处理起来最准确
        return json.dumps(table_names, ensure_ascii=False)

    except Exception as e:
        return f"Error listing tables: {str(e)}"


@mcp.tool(exclude_args=["plugin_config"])
def get_database_schema(
    table_names: List[str] | None = None,
    plugin_config: dict[str, str] = {},
) -> str:
    """
    获取数据库 Schema。

    Args:
        table_names: 可选。指定需要查询的表名列表。如果不填，默认返回所有表名（不含列详情），以节省 token。
    """
    try:
        engine = manager.get_engine(plugin_config)
        inspector = inspect(engine)
        all_tables = inspector.get_table_names()

        # 策略：如果不指定表名，只返回表名列表（防止 Token 爆炸）
        if not table_names:
            return f"Found tables: {', '.join(all_tables)}\n\nPlease call this tool again with specific `table_names` to get column details."

        # 获取指定表的详细信息
        schema_info = []
        for table in table_names:
            if table not in all_tables:
                schema_info.append(f"Table '{table}' not found.")
                continue

            columns = inspector.get_columns(table)
            # 格式化为更易读的 Prompt 格式
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
def execute_sql(query: str, plugin_config: dict[str, str] = {}) -> str | dict:
    """
    Args:
        query: 要执行的 SQL Select 语句。
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
            keys = result.keys()
            rows = result.fetchmany(50)  # 限制返回行数

            if not rows:
                return "Query returned no results."

            output = [dict(zip(keys, row)) for row in rows]
            return {"model_text": json.dumps(output, default=str, ensure_ascii=False)}

    except Exception as e:
        return f"SQL Execution Error: {str(e)}"


if __name__ == "__main__":
    mcp.run()
