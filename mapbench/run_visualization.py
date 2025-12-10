import scanpy as sc
import os
import matplotlib.pyplot as plt
from mapbench.metrics import uncertainty, DALogFC, MapQC
import mapqc
import milopy.utils
import milopy.plot as milopl

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection

pd.DataFrame.iteritems = pd.DataFrame.items


def plot_uncertainty(adata, output_dir="./"):
    """Plot cell uncertainty scores with check for pre-computed values"""
    if 'uncertainty' not in adata.uns:
        print("→ Computing uncertainty...")
        adata = uncertainty(adata)
    else:
        print("→ Using pre-computed uncertainty from adata.uns")
    print("→ Plotting uncertainty...")
    adata_query = adata[adata.obs["ref_query"] == "query"].copy()
    adata_query.obs["uncertainty"] = adata.uns['uncertainty']
    # sc.pl.umap 绘图，但 fig=plt.gcf() 
    sc.pl.umap(adata_query, color="uncertainty", title="Uncertainty", show=False, frameon=False, cmap='RdBu_r')
    
    # 获取当前图形对象
    fig = plt.gcf()
    
    # 1. 栅格化点位以减小 PDF 文件大小
    for ax in fig.get_axes():
        for child in ax.get_children():
            # 识别散点图的点集 (PathCollection)
            if isinstance(child, PathCollection):
                child.set_rasterized(True) # 启用栅格化

    # 2. 保存时使用 bbox_inches='tight' 修正图例截断并设置 DPI
    plt.savefig(
        os.path.join(output_dir, "uncertainty_plot.pdf"),
        bbox_inches='tight', # 确保图例不被截断
        dpi=300              # 设置栅格化图像的分辨率
    )
    print("✅ Saved: uncertainty_plot.pdf")


def plot_DAlogFC(adata, output_dir=".", **kwargs):
    """Plot cell DAlogFC scores with check for pre-computed values"""
    if 'nhood_adata' not in adata.uns:
        print("→ Computing DAlogFC...")
        adata = DALogFC(adata)
    else:
        print("→ Using pre-computed DAlogFC from adata.uns")

    print("→ Building neighborhood graph...")
    milopy.utils.build_nhood_graph(adata)

    print("→ Plotting DAlogFC...")
    plt.rcParams["figure.figsize"] = [10, 10]
    milopl.plot_nhood_graph(
        adata,
        show=False,
        min_size=1,
        **kwargs
    )
    
    # 获取当前图形对象
    fig = plt.gcf()
    
    # 1. 栅格化点位以减小 PDF 文件大小
    # 尽管 milopl.plot_nhood_graph 主要绘制邻域，但为了覆盖可能的背景点或确保其他 PathCollection 元素被处理
    for ax in fig.get_axes():
        for child in ax.get_children():
            # 识别散点图的点集 (PathCollection)
            if isinstance(child, PathCollection):
                child.set_rasterized(True) # 启用栅格化

    # 2. 保存时使用 bbox_inches='tight' 修正图例截断并设置 DPI
    plt.savefig(
        os.path.join(output_dir, "DAlogFC_plot.pdf"),
        bbox_inches='tight', # 确保图例不被截断
        dpi=300              # 设置栅格化图像的分辨率
    )
    print("✅ Saved: DAlogFC_plot.pdf")



def plot_mapQC(adata, output_dir="."):
    if 'mapqc_score' not in adata.obs:
        print("→ mapQC scores not found in adata.obs, computing...")
        _ = MapQC(
            adata,
            embedding="X",
            ref_query_key="ref_query",
            ref_key="ref",
            query_key="query",
            study_key="dataset_id",
            n_nhoods=300,
            k_min=1500,
            k_max=10000,
            seed=10,
            overwrite=True,
            return_nhood_info_df=False,
            return_sample_dists_to_ref_df=False,
        )
        # 添加 case/control 标签（mapQC 用）
        target_celltype = adata.uns["oor_celltype"]
        adata.obs['case_control'] = None
        adata.obs.loc[adata.obs['ref_query'] == 'query', 'case_control'] = 'control'
        adata.obs.loc[
            (adata.obs['ref_query'] == 'query') & (adata.obs['cell_type'] == target_celltype),
            'case_control'
        ] = 'case'
    else:
        print("→ Using pre-computed mapQC scores from adata.obs")
    # _ = mapqc.evaluate(
    #     adata,
    #     case_control_key="case_control",
    #     case_cats=["case"],
    #     control_cats=["control"]
    # )

    print("→ Plotting mapQC...")
    fig = mapqc.pl.umap.mapqc_scores_binary(adata, return_fig=True)

    # 1. 栅格化点位以减小 PDF 文件大小
    for ax in fig.get_axes():
        for child in ax.get_children():
            # 识别散点图的点集
            if isinstance(child, PathCollection):
                child.set_rasterized(True) # 启用栅格化

    # 2. 保存时使用 bbox_inches='tight' 修正图例截断
    fig.savefig(
        os.path.join(output_dir, "mapQC_plot.pdf"), 
        bbox_inches='tight', # 确保图例不被截断
        dpi=300              # 设置栅格化图像的分辨率
    ) 

    print("✅ Saved: mapQC_plot.pdf")



# def plot_all_methods(adata_path, output_dir="./"):
#     """
#     执行 uncertainty / DAlogFC / mapQC 的计算 + 可视化。
#     每个函数内部自带计算逻辑。
#     """
#     adata = sc.read_h5ad(adata_path)
#     adata_name = os.path.basename(adata_path).replace(".h5ad", "")
#     if "oor_celltype" not in adata.uns:
#         raise ValueError("`adata.uns['oor_celltype']` is required.")
#     if "ref_query" not in adata.obs:
#         raise ValueError("`adata.obs['ref_query']` is required.")


#     print(f"=== Visualizing {adata_name} ===")
#     plot_uncertainty(adata, output_dir)
#     plot_DAlogFC(adata, output_dir)
#     plot_mapQC(adata, output_dir)
#     print(f"✅ All visualizations complete for {adata_name}.")


def plot_all_methods(adata, output_dir="./"):
    """
    执行 uncertainty / DAlogFC / mapQC 的计算 + 可视化。
    每个函数内部自带计算逻辑。
    """
    if "oor_celltype" not in adata.uns:
        raise ValueError("`adata.uns['oor_celltype']` is required.")
    if "ref_query" not in adata.obs:
        raise ValueError("`adata.obs['ref_query']` is required.")


    print(f"=== Visualizing  ===")
    plot_uncertainty(adata, output_dir)
    plot_DAlogFC(adata, output_dir)
    plot_mapQC(adata, output_dir)
    print(f"✅ All visualizations complete .")



def plot_auprc_summary(csv_path="evaluation_result.csv", output_dir="./"):
    """
    从 evaluation_result.csv 读取 AUPRC 数据，绘制 boxplot/stripplot。
    每一行是一个样本，每一列是一个方法。
    """
    df = pd.read_csv(csv_path, sep='\t', header=0)

    if not all(col in df.columns for col in ["Uncertainty", "DAlogFC", "mapQC"]):
        raise ValueError("Missing required AUPRC columns in CSV.")

    # melt成长表格：method、AUPRC、adata_name
    df_melted = df.melt(id_vars=["adata_name"], 
                        value_vars=["Uncertainty", "DAlogFC", "mapQC"],
                        var_name="Method", value_name="AUPRC")

    plt.figure(figsize=(6, 4))
    sns.boxplot(data=df_melted, x="Method", y="AUPRC", palette="Set2", width=0.6, fliersize=0)
    sns.stripplot(data=df_melted, x="Method", y="AUPRC", color="black", size=4, jitter=True)

    plt.title("AUPRC Summary across Methods", fontsize=13)
    plt.ylabel("AUPRC", fontsize=11)
    plt.xlabel("")
    plt.ylim(0, 1.05)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "AUPRC_summary_plot.pdf"))
    print("✅ Saved: AUPRC_summary_plot.pdf")


