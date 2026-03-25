# dingent-hub

A curated marketplace of [FastMCP](https://github.com/jlowin/fastmcp) plugins that extend AI assistants with domain-specific tools — from bioinformatics and single-cell RNA-seq analysis to natural-language SQL and retrieval-augmented generation (RAG).

---

## Table of Contents

- [Overview](#overview)
- [Plugins](#plugins)
  - [Bio Data Loader](#bio-data-loader)
  - [Ding RAG](#ding-rag)
  - [Ding Text2SQL](#ding-text2sql)
  - [Single Cell Analyzer](#single-cell-analyzer)
- [Tech Stack](#tech-stack)
- [Getting Started](#getting-started)
- [Plugin Configuration](#plugin-configuration)
- [Contributing](#contributing)
- [License](#license)

---

## Overview

`dingent-hub` is a plugin hub built on the [Model Context Protocol (MCP)](https://modelcontextprotocol.io/). Each plugin is an independent MCP server that exposes a set of typed tools to any compatible AI client. Plugins are self-contained Python projects managed with [uv](https://github.com/astral-sh/uv) and run via stdio transport.

```
dingent-hub/
├── market.json               # Hub metadata (version, plugin counts)
├── LICENSE
└── plugins/
    ├── bio-data-loader/      # Download & convert bioinformatics data
    ├── ding-rag/             # RAG with Chroma vector database
    ├── ding-text2sql/        # Natural language → SQL
    └── single-cell-analyzer/ # Single-cell RNA-seq analysis pipeline
```

---

## Plugins

### Bio Data Loader

> **Version:** 1.0.0 · **Python:** ≥ 3.14

Download and convert bioinformatics data files into analysis-ready formats.

| Tool | Description |
|------|-------------|
| `fetch_file_from_url` | Stream-download a file from a URL; auto-extracts ZIP / TAR.GZ archives |
| `convert_10x_to_h5ad` | Convert 10x Genomics output (matrix.mtx + barcodes + features) to H5AD |
| `convert_csv_to_h5ad` | Convert a CSV expression matrix to H5AD (with optional transpose) |
| `load_demo_dataset` | Load built-in demo datasets (`pbmc3k`, `paul15`) for quick testing |

**Key dependencies:** `fastmcp`, `scanpy`, `ipython`

---

### Ding RAG

> **Version:** 1.1.0 · **Python:** ≥ 3.12

Provide retrieval-augmented generation using a [Chroma](https://www.trychroma.com/) vector database.

| Tool | Description |
|------|-------------|
| `search_knowledge` | Query the knowledge base and return the most relevant document chunks |

**Supported embedding providers:** OpenAI, HuggingFace, Google, Cohere, SentenceTransformer, Chroma default

**Key dependencies:** `fastmcp`, `langchain-chroma`, `langchain-community`, `openai`

---

### Ding Text2SQL

> **Version:** 1.3.0 · **Python:** ≥ 3.10

Convert natural-language questions into SQL queries and execute them against a relational database.

| Tool | Description |
|------|-------------|
| `list_tables` | List all tables available in the connected database |
| `get_database_schema` | Return a compact, LLM-optimised schema representation |
| `execute_sql` | Execute a `SELECT` query (destructive statements are blocked) |

**Security:** Only `SELECT` queries are allowed. `DROP`, `DELETE`, `UPDATE`, `INSERT`, `ALTER`, `GRANT`, and `TRUNCATE` are all rejected. SQLite URIs are disabled.

**Key dependencies:** `fastmcp`, `sqlmodel`, `pymysql`

---

### Single Cell Analyzer

> **Version:** 1.0.0 · **Python:** ≥ 3.14

End-to-end single-cell RNA-seq analysis pipeline with built-in visualisation.

| Tool | Description |
|------|-------------|
| `quality_control_analysis` | Compute QC metrics (genes/cell, UMI counts, mitochondrial %) and filter cells |
| `run_clustering_and_umap` | Normalise → log-transform → HVG selection → PCA → Leiden clustering → UMAP |
| `find_marker_genes` | Identify differentially expressed genes per cluster (Wilcoxon); generates dot plots |
| `run_paga_trajectory` | PAGA trajectory analysis with ForceAtlas2 embedding |
| `generate_markdown_report` | Produce a Markdown report with embedded Base64 images |

**Key dependencies:** `fastmcp`, `scanpy`, `anndata`, `pandas`, `matplotlib`, `seaborn`, `leidenalg`

---

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Protocol | [Model Context Protocol (MCP)](https://modelcontextprotocol.io/) |
| Framework | [FastMCP](https://github.com/jlowin/fastmcp) |
| Language | Python 3.10 – 3.14 |
| Package manager | [uv](https://github.com/astral-sh/uv) |
| Transport | stdio |
| Vector DB (RAG) | [Chroma](https://www.trychroma.com/) |
| Bioinformatics | [Scanpy](https://scanpy.readthedocs.io/) / [AnnData](https://anndata.readthedocs.io/) |
| Database ORM | [SQLModel](https://sqlmodel.tiangolo.com/) / SQLAlchemy |

---

## Getting Started

### Prerequisites

- **Python 3.10+** (exact version depends on the plugin — see above)
- **[uv](https://docs.astral.sh/uv/getting-started/installation/)** package manager

### Running a plugin

Each plugin is a standalone MCP server. Navigate into the plugin directory and start it with `uv`:

```bash
cd plugins/ding-text2sql
uv run main.py --verbose
```

`uv` automatically creates an isolated virtual environment and installs all locked dependencies from `uv.lock` on first run.

### Connecting to an MCP client

Add the plugin to your MCP client configuration. Example for `ding-text2sql`:

```json
{
  "mcpServers": {
    "ding-text2sql": {
      "command": "uv",
      "args": ["run", "main.py", "--verbose"],
      "cwd": "/path/to/plugins/ding-text2sql",
      "env": {
        "DB_URI": "mysql+pymysql://user:password@host/dbname"
      }
    }
  }
}
```

---

## Plugin Configuration

Every plugin is configured via a `plugin.toml` file and exposes a `config_schema` that the MCP client passes as environment variables or structured config at startup.

### ding-text2sql

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `db_uri` | string (secret) | ✅ | SQLAlchemy database URL — see [docs](https://docs.sqlalchemy.org/en/20/core/engines.html#database-urls) |

### ding-rag

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `client_type` | string | ✅ | Chroma client mode: `persistent`, `http`, or `ephemeral` |
| `collection_name` | string | ✅ | Chroma collection to query |
| `embedding_provider` | string | ✅ | Embedding backend (e.g. `openai`, `huggingface`) |
| `embedding_model_name` | string | ✅ | Model name for the chosen embedding provider |
| `client_path` | string | — | Local path for `persistent` client |
| `client_host` | string | — | Host for `http` client |
| `client_port` | integer | — | Port for `http` client |
| `embedding_api_key` | string (secret) | — | API key for the embedding provider |
| `embedding_api_base` | string (secret) | — | Custom API base URL |
| `embedding_extra_kwargs` | dict | — | Additional provider-specific kwargs |
| `default_n_results` | integer | — | Number of chunks to return (default: 3) |

---

## Contributing

Contributions are welcome! To add a new plugin:

1. Create a new directory under `plugins/` using the existing plugins as a template.
2. Implement your MCP server in `main.py` using FastMCP.
3. Declare metadata and configuration schema in `plugin.toml`.
4. Pin dependencies in `pyproject.toml` and generate `uv.lock` with `uv lock`.
5. Update `market.json` to reflect the new plugin count.
6. Open a pull request with a description of your plugin's purpose and tools.

---

## License

This project is licensed under the [MIT License](LICENSE).
