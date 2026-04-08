# biomarker-MCP

Natural language interface for celltype marker query through MCP.

## 🪩 What can it do?

- query celltype marker from CellMarker database by natural language
- will add more database

## ❓ Who is this for?

- Anyone who wants to do query celltype marker with natural language!
- Agent developers who want to query cell markers for their applications

## 🌐 Where to use it?

You can use biomarker-mcp in most AI clients, plugins, or agent frameworks that support the MCP:

- AI clients, like Cherry Studio
- Plugins, like Cline
- Agent frameworks, like Agno 

## 🎬 Demo

A demo showing query celltype markers in a AI client Cherry Studio using natural language based on biomarker-mcp

https://github.com/user-attachments/assets/71268f6f-c74d-4142-ad7a-893b411d748a


## 📚 Documentation

scmcphub's complete documentation is available at https://docs.scmcphub.org

## 🏎️ Quickstart

### Install

Install from PyPI
```
pip install biomarker-mcp
```
you can test it by running
```
biomarker-mcp run
```


#### run biomarker-mcp locally
Refer to the following configuration in your MCP client:

check path
```
$ which biomarker 
/home/test/bin/biomarker-mcp
```

```
"mcpServers": {
  "biomarker-mcp": {
    "command": "/home/test/bin/biomarker-mcp",
    "args": [
      "run"
    ]
  }
}
```

#### run biomarker-server remotely
Refer to the following configuration in your MCP client:

run it in your server
```
biomarker-mcp run --transport shttp --port 8000
```

Then configure your MCP client in local AI client, like this:
```

"mcpServers": {
  "biomarker-mcp": {
    "url": "http://localhost:8000/mcp"
  }
}
```

## 🤝 Contributing

If you have any questions, welcome to submit an issue, or contact me(hsh-me@outlook.com). Contributions to the code are also welcome!

## Citing

If you use biomarker-mcp in for your research, please consider citing  following work: 
> Congxue Hu, Tengyue Li, Yingqi Xu, Xinxin Zhang, Feng Li, Jing Bai, Jing Chen, Wenqi Jiang, Kaiyue Yang, Qi Ou, Xia Li, Peng Wang, Yunpeng Zhang, CellMarker 2.0: an updated database of manually curated cell markers in human/mouse and web tools based on scRNA-seq data, Nucleic Acids Research, Volume 51, Issue D1, 6 January 2023, Pages D870–D876, https://doi.org/10.1093/nar/gkac947
