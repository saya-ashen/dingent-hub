from typing import Any
from fastmcp import FastMCP
import scanpy as sc
import matplotlib.pyplot as plt
import io
import os
import base64
import json

# 创建 MCP Server
mcp = FastMCP("Single-Cell-Analyzer")

# --- 1. 全局图片缓存 (关键修改) ---
# 用于在内存中临时存储 Base64 图片数据，避免 LLM 上下文溢出
IMAGE_STORE = {}

# 全局绘图设置
sc.settings.set_figure_params(
    dpi=100, frameon=False, vector_friendly=True, color_map="viridis"
)
plt.switch_backend("Agg")


def get_top_markers_text(adata, n_top=5):
    """提取每个簇的前 n 个 Marker 基因，返回文本供 LLM 阅读"""
    result = {}
    groups = adata.uns["rank_genes_groups"]["names"].dtype.names
    for group in groups:
        genes = [
            str(adata.uns["rank_genes_groups"]["names"][i][group]) for i in range(n_top)
        ]
        result[group] = genes
    return json.dumps(result, indent=2)


def save_plot_to_store(image_key: str) -> str:
    """
    将当前 Matplotlib 图片转换为 Base64，存入全局缓存，并返回 Base64 字符串。
    Args:
        image_key: 图片的唯一标识符 (例如 'qc_plot')
    """
    buf = io.BytesIO()
    plt.savefig(buf, format="png", bbox_inches="tight")
    buf.seek(0)
    b64_str = base64.b64encode(buf.read()).decode("utf-8")
    plt.close()

    # 构造完整的 Data URI
    full_b64 = f"data:image/png;base64,{b64_str}"

    # 存入全局缓存
    IMAGE_STORE[image_key] = full_b64

    return full_b64


def format_response(
    model_text: str, display_title: str, markdown_content: str
) -> dict[str, Any]:
    return {
        "model_text": model_text,
        "display": [
            {"type": "markdown", "title": display_title, "content": markdown_content}
        ],
    }


@mcp.tool()
def quality_control_analysis(file_path: str) -> dict[str, Any]:
    """
    执行质量控制 (QC)。
    分析完成后，图片会被缓存为 ID: 'qc_plot'。
    """
    if not os.path.exists(file_path):
        return {"model_text": "Error: File not found.", "display": []}

    adata = sc.read_h5ad(file_path)

    # 计算指标
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # 绘图
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )

    # --- 保存到缓存，Key 为 'qc_plot' ---
    img_b64 = save_plot_to_store("qc_plot")

    # 数据过滤
    original_cells = adata.n_obs
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]
    remaining_cells = adata.n_obs

    new_path = file_path.replace(".h5ad", "_qc.h5ad")
    adata.write(new_path)

    # --- 告诉 LLM 图片 ID ---
    model_msg = (
        f"QC Analysis completed.\n"
        f"Original: {original_cells}, Remaining: {remaining_cells}.\n"
        f"Filtered data saved to: {new_path}.\n"
        f"IMPORTANT: The QC plot has been cached with ID 'qc_plot'. "
        f"When generating the report, use the placeholder {{{{qc_plot}}}} to insert it."
    )

    display_content = (
        f"### 🧬 质量控制分析结果\n"
        f"- **过滤前**: {original_cells}\n"
        f"- **过滤后**: {remaining_cells}\n"
        f"![QC Plot]({img_b64})"
    )

    return format_response(model_msg, "QC Analysis Result", display_content)


@mcp.tool()
def run_clustering_and_umap(file_path: str) -> dict[str, Any]:
    """
    运行聚类。
    分析完成后，图片会被缓存为 ID: 'umap_plot'。
    """
    adata = sc.read_h5ad(file_path)

    # (标准分析流程，简化展示)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    # 绘图
    sc.pl.umap(adata, color=["leiden"], title="Cell Clusters", show=False)

    # --- 保存到缓存，Key 为 'umap_plot' ---
    img_b64 = save_plot_to_store("umap_plot")

    num_clusters = len(adata.obs["leiden"].unique())
    new_path = file_path.replace(".h5ad", "_processed.h5ad")
    adata.write(new_path)

    model_msg = (
        f"Clustering completed. Found {num_clusters} clusters.\n"
        f"IMPORTANT: The UMAP plot has been cached with ID 'umap_plot'. "
        f"When generating the report, use the placeholder {{{{umap_plot}}}} to insert it."
        f" Data saved to: {new_path}."
    )

    display_content = (
        f"### 🗺️ 聚类结果 (UMAP)\n"
        f"共发现 **{num_clusters}** 个细胞簇。\n"
        f"![UMAP Plot]({img_b64})"
    )

    return format_response(model_msg, "Clustering Visualization", display_content)


@mcp.tool()
def find_marker_genes(file_path: str, groupby: str = "leiden") -> dict[str, Any]:
    """
    计算差异表达基因 (Marker Genes)，用于鉴定细胞类型。
    分析完成后，Dotplot 图片会被缓存为 ID: 'marker_plot'。
    """
    adata = sc.read_h5ad(file_path)

    # 确保已经做过聚类
    if groupby not in adata.obs:
        return {
            "model_text": f"Error: '{groupby}' not found. Run clustering first.",
            "display": [],
        }

    # 计算差异基因 (Wilcoxon rank-sum)
    sc.tl.rank_genes_groups(adata, groupby, method="wilcoxon")

    # 绘图：Dotplot 是最直观的 Marker 展示方式
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, show=False)
    img_b64 = save_plot_to_store("marker_plot")

    # 提取 Top 基因列表给 LLM
    top_genes_json = get_top_markers_text(adata, n_top=5)

    # 保存结果
    new_path = file_path.replace(".h5ad", "_markers.h5ad")
    adata.write(new_path)

    model_msg = (
        f"Marker gene analysis completed.\n"
        f"Top 5 genes per cluster identified: {top_genes_json}\n"  # 把基因直接给 LLM，让 LLM 进行生物学注释
        f"IMPORTANT: The Dotplot has been cached as 'marker_plot'. Use {{marker_plot}} in the report.\n"
        f"TASK FOR LLM: Based on the gene list above, please infer the cell type for each cluster in your response."
        f" Data saved to: {new_path}."
    )

    display_content = (
        f"### 🧬 差异基因分析 (Markers)\n"
        f"已通过 Wilcoxon 检验计算各簇特征基因。\n"
        f"![Marker Dotplot]({img_b64})\n\n"
        f"**Top Genes per Cluster:**\n```json\n{top_genes_json}\n```"
    )

    return format_response(model_msg, "Marker Analysis Result", display_content)


@mcp.tool()
def run_trajectory_analysis(file_path: str) -> dict[str, Any]:
    """
    运行 PAGA (Partition-based Graph Abstraction) 进行细胞轨迹推断。
    适用于有发育关系的数据集。缓存图片 ID: 'paga_plot'。
    """
    adata = sc.read_h5ad(file_path)

    # 确保已有 Neighbors 和 Leiden 结果
    if "leiden" not in adata.obs:
        return {"model_text": "Error: Run clustering first.", "display": []}

    # 运行 PAGA
    sc.tl.paga(adata, groups="leiden")

    # 绘图：PAGA 拓扑图 + UMAP 嵌入
    sc.pl.paga(adata, show=False)
    # 也可以结合 UMAP 绘制 (sc.tl.draw_graph 比较慢，这里只画拓扑结构)

    img_b64 = save_plot_to_store("paga_plot")

    new_path = file_path.replace(".h5ad", "_paga.h5ad")
    adata.write(new_path)

    model_msg = (
        f"PAGA trajectory analysis completed.\n"
        f"Connectivity graph generated. Connectivity threshold indicates the strength of relation between clusters.\n"
        f"Cached image ID: 'paga_plot'. Use {{paga_plot}} in report."
        f" Data saved to: {new_path}."
    )

    display_content = (
        f"### 🕸️ 细胞轨迹推断 (PAGA)\n"
        f"展示了各细胞簇之间的拓扑连接关系（线条越粗表示连通性越强）。\n"
        f"![PAGA Plot]({img_b64})"
    )

    return format_response(model_msg, "Trajectory Analysis", display_content)


@mcp.tool()
def generate_markdown_report(report_title: str, markdown_body: str) -> dict[str, Any]:
    """
    生成包含 Base64 图片的 Markdown 报告。

    Args:
        report_title: 报告标题
        markdown_body: 报告正文。
                       **关键**: 如果需要插入图片，请在文本中使用 {{image_id}} 占位符。
                       例如: "这是 QC 结果: {{qc_plot}}" 或 "这是聚类图: {{umap_plot}}"。
                       工具会自动将其替换为 Base64 图片代码。
    """

    # 1. 替换占位符
    # 我们遍历缓存中的图片，查找 markdown_body 中是否有对应的占位符 {{key}}
    # 如果有，替换为标准的 Markdown 图片语法 ![key](base64_data)

    processed_body = markdown_body

    for key, b64_data in IMAGE_STORE.items():
        placeholder = f"{{{{{key}}}}}"  # 匹配字符串 "{{key}}"
        if placeholder in processed_body:
            # 替换为 Markdown 图片语法
            markdown_image = f"![{key}]({b64_data})"
            processed_body = processed_body.replace(placeholder, markdown_image)

    # 2. 组装最终 Markdown 内容
    final_content = f"# {report_title}\n\n{processed_body}\n\n---\n*Generated by Single-Cell-Analyzer MCP*"

    # 3. 保存到本地文件 (Base64 很大，建议保存为文件查看)

    # 4. 构建返回信息
    # 注意：我们不在 model_text 里返回整个 Base64 内容，防止刷屏。
    model_msg = (
        f"Report has been generated. All placeholders replaced with Base64 images."
    )

    # 在前端展示部分，我们可以展示一个缩略版本，或者直接提示文件已生成

    return format_response(model_msg, "Report Generated", final_content)


if __name__ == "__main__":
    mcp.run()
