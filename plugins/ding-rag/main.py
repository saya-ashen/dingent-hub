from typing import Any, Literal
import chromadb
from chromadb.utils import embedding_functions
from fastmcp import FastMCP
from pydantic import BaseModel

# 定义 Client 类型
ClientType = Literal["persistent", "http", "ephemeral"]


class PluginConfig(BaseModel):
    client_type: ClientType = "persistent"
    client_path: str = "./genbase_chroma_db/"  # 仅 persistent 模式使用
    client_host: str = "localhost"  # 仅 http 模式使用
    client_port: int = 8000  # 仅 http 模式使用

    # --- Collection & Embedding 配置 ---
    collection_name: str
    embedding_provider: str
    embedding_model_name: str | None = None
    embedding_api_key: str | None = None
    embedding_api_base: str | None = None
    embedding_extra_kwargs: dict[str, Any] | None = None
    default_n_results: int = 3


mcp = FastMCP("rag", stateless_http=True)

# 全局缓存池
# Key: client_config_hash -> ClientObj
_clients_cache: dict[str, chromadb.ClientAPI] = {}

# Key: client_hash|collection_name|provider|model_name -> CollectionObj
_collections_cache: dict[str, chromadb.Collection] = {}


def get_ef(provider: str, **kwargs):
    """
    (保持不变) 根据 provider 名称动态返回 Chroma 内置的 EmbeddingFunction
    """
    provider_map = {
        "openai": embedding_functions.OpenAIEmbeddingFunction,
        "huggingface": embedding_functions.HuggingFaceEmbeddingFunction,
        "google": embedding_functions.GooglePalmEmbeddingFunction,
        "cohere": embedding_functions.CohereEmbeddingFunction,
        "sentence_transformer": embedding_functions.SentenceTransformerEmbeddingFunction,
        "default": embedding_functions.DefaultEmbeddingFunction,
    }

    ef_class = provider_map.get(provider.lower())
    if not ef_class:
        raise ValueError(f"不支持的 provider: {provider}")

    valid_kwargs = {k: v for k, v in kwargs.items() if v is not None}

    try:
        return ef_class(**valid_kwargs)
    except TypeError as e:
        if "model_name" in valid_kwargs and "unexpected keyword argument" in str(e):
            return ef_class(model_name=valid_kwargs["model_name"])
        raise e


def get_client_safe(config: PluginConfig) -> tuple[str, chromadb.ClientAPI]:
    """
    根据配置获取 Client，并返回 (client_hash, client_instance)
    """
    # 构造 Client 的唯一标识 Key
    if config.client_type == "persistent":
        cache_key = f"persistent|{config.client_path}"
    elif config.client_type == "http":
        cache_key = f"http|{config.client_host}:{config.client_port}"
    else:
        cache_key = "ephemeral"  # 内存模式通常不复用或单例，视需求而定

    # 如果缓存中有，直接返回
    if cache_key in _clients_cache:
        return cache_key, _clients_cache[cache_key]

    # 初始化 Client
    print(f"正在初始化 Chroma Client: {cache_key}...")
    if config.client_type == "persistent":
        client = chromadb.PersistentClient(path=config.client_path)
    elif config.client_type == "http":
        client = chromadb.HttpClient(
            host=config.client_host, port=str(config.client_port)
        )
    else:
        client = chromadb.EphemeralClient()

    _clients_cache[cache_key] = client
    return cache_key, client


def get_collection_safe(
    client: chromadb.ClientAPI, client_key: str, name: str, provider: str, **ef_kwargs
):
    """
    获取 Collection，Key 必须包含 client_key 以区分不同数据源
    """
    model_name = ef_kwargs.get("model_name", "default")

    cache_key = f"{client_key}|{name}|{provider}|{model_name}"

    if cache_key in _collections_cache:
        return _collections_cache[cache_key]

    # 实例化 EF
    embedding_fn = get_ef(provider, **ef_kwargs)

    try:
        col = client.get_collection(name=name, embedding_function=embedding_fn)
    except Exception:
        # 如果集合不存在，返回 None，交由上层处理
        return None

    _collections_cache[cache_key] = col
    return col


@mcp.tool(exclude_args=["plugin_config"])
def search_knowledge(
    query: str,
    n_results: int | str | None = 3,
    plugin_config: PluginConfig | None = None,
) -> dict[str, Any]:
    """
    从本地知识库中检索相关信息。支持动态配置数据库路径。
    """
    assert plugin_config, "PluginConfig is required"

    # 1. 解包配置
    n_results = int(n_results or plugin_config.default_n_results)

    ef_params = {
        "api_key": plugin_config.embedding_api_key,
        "model_name": plugin_config.embedding_model_name,
        "api_base": plugin_config.embedding_api_base,
    }
    ef_params.update(plugin_config.embedding_extra_kwargs or {})

    try:
        # 2. 获取动态 Client
        client_key, client = get_client_safe(plugin_config)

        # 3. 获取集合 (传入 client 和 client_key)
        collection = get_collection_safe(
            client=client,
            client_key=client_key,
            name=plugin_config.collection_name,
            provider=plugin_config.embedding_provider,
            **ef_params,
        )

        if collection is None:
            return {
                "model_text": f"错误: 在数据库 '{client_key}' 中未找到知识库 '{plugin_config.collection_name}'。",
                "display": {},
            }

        # 4. 执行查询
        results = collection.query(query_texts=[query], n_results=n_results)

        documents = results["documents"][0] if results["documents"] else []
        metadatas = results["metadatas"][0] if results["metadatas"] else []

        if not documents:
            return {"model_text": "未找到相关信息。", "display": {}}

        response_text = f"在 '{plugin_config.collection_name}' 中找到以下信息:\n\n"
        for i, doc in enumerate(documents):
            meta = metadatas[i] if i < len(metadatas) else {}
            # 处理 metadata 为 None 的情况
            meta = meta or {}
            title = meta.get("title", "Unknown Source")
            response_text += f"--- 来源: {title} ---\n{doc}\n\n"

        return {"model_text": response_text.strip(), "display": {}}

    except Exception as e:
        import traceback

        traceback.print_exc()
        return {"model_text": f"检索过程发生系统错误: {str(e)}", "display": {}}


if __name__ == "__main__":
    mcp.run()
