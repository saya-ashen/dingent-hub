# ChEMBL-MCP-Server
[![smithery badge](https://smithery.ai/badge/@JackKuo666/chembl-mcp-server)](https://smithery.ai/server/@JackKuo666/chembl-mcp-server)

A FastMCP wrapper server based on the chembl_webresource_client package, providing API access to the ChEMBL database.

## Features

- Complete API access to the ChEMBL database
- Asynchronous API calls implemented using FastMCP framework
- Built-in error handling and timeout mechanisms
- Support for both HTTP and stdio transport methods
- Complete type annotations and docstrings

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/ChEMBL-MCP-Server.git
cd ChEMBL-MCP-Server

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Starting the Server

```bash
# Start HTTP server with default configuration
python chembl_searver.py

# Specify host and port
python chembl_searver.py --host 0.0.0.0 --port 8080

# Use stdio transport
python chembl_searver.py --transport stdio

# Set log level
python chembl_searver.py --log-level DEBUG
```

### Available Parameters

- `--host`: Server host address, defaults to 127.0.0.1
- `--port`: Server port, defaults to 8000
- `--transport`: Transport method, choose between http or stdio, defaults to http
- `--log-level`: Log level, choose from DEBUG, INFO, WARNING, ERROR, CRITICAL, defaults to INFO

## API Functions

The server provides the following API functions:

### Data Entity APIs

- `example_activity`: Get activity data
- `example_assay`: Get assay data
- `example_target`: Get target data
- `example_molecule`: Get molecule data
- `example_drug`: Get drug data
- More data entity APIs...

### Chemical Tool APIs

- `example_canonicalizeSmiles`: Canonicalize SMILES strings
- `example_smiles2inchi`: Convert SMILES to InChI
- `example_smiles2svg`: Convert SMILES to SVG image
- `example_structuralAlerts`: Get structural alerts
- More chemical tool APIs...

## Examples

Check the `chembl_search.py` file for examples of using various APIs.

## Dependencies

- chembl_webresource_client: ChEMBL Web Service Client
- mcp: MCP Framework
- fastapi: FastAPI Framework
- uvicorn: ASGI Server
- asyncio: Asynchronous I/O Library

## License

[MIT](LICENSE)