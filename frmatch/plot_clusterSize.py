
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib

def plot_clusterSize(adata_E1, cluster_header_E1, proportion = False, adata_E2 = None, cluster_header_E2 = None, name_E1 = "E1", name_E2 = "E2", 
                     color_E1 = "#F0E442", color_E2 = "#56B4E9", y_offset = 0.05, width = 10, height = 10, save = False): 
    """\
    Plots cluster size for input AnnData objects.

    Parameters
    ----------
        adata_E1: AnnData
            Annotated data matrix.
        cluster_header_E1: str
            Column in `adata_E1.obs` storing cell annotation.
        adata_E2: AnnData (default: None)
            Annotated data matrix.
        cluster_header_E2: str (default: None)
            Column in `adata_E2.obs` storing cell annotation.
        name_E1: str (default: "E1")
            Label for `adata_E1`.
        name_E2: str (default: "E2")
            Label for `adata_E2`.
        color_E1: str (default: yellow)
            Bar color for `adata_E1`.
        color_E2: str (default: blue)
            Bar color for `adata_E2`.
        width: int (default: 10)
            Width of plot.
        height: int (default: 10)
            Height of plot.
    """
    if adata_E2 and cluster_header_E2: 
        fig, axes = plt.subplots(2, figsize = (width, height))
        title = f"{name_E1} ({len(np.unique(adata_E1.obs[cluster_header_E1]))} clusters, {adata_E1.shape[0]} cells)"
        dictionary = dict(adata_E1.obs[cluster_header_E1].value_counts())
        if proportion: 
            total = sum(dictionary.values())
            for key in dictionary: dictionary[key] = round(dictionary[key]/total, 3)
        p = axes[0].bar(dictionary.keys(), dictionary.values(), color = color_E1)
        ax.set_ylim(0, list(dictionary.values())[0] + y_offset)
        axes[0].bar_label(p, fontsize = 8)
        axes[0].set_title(title)
        if proportion: axes[0].set_ylabel("Proportion")
        else: axes[0].set_ylabel("Size")
        axes[0].set_xlabel("Cluster")
        axes[0].set_xticklabels(dictionary.keys(), rotation = 90)

        title = f"{name_E2} ({len(np.unique(adata_E2.obs[cluster_header_E2]))} clusters, {adata_E2.shape[0]} cells)"
        dictionary = dict(adata_E2.obs[cluster_header_E2].value_counts())
        if proportion: 
            total = sum(dictionary.values())
            for key in dictionary: dictionary[key] = round(dictionary[key]/total, 3)
        p = axes[1].bar(dictionary.keys(), dictionary.values(), color = color_E2)
        axes[1].bar_label(p, fontsize = 8)
        axes[1].set_title(title)
        if proportion: axes[1].set_ylabel("Proportion")
        else: axes[1].set_ylabel("Size")
        axes[1].set_xlabel("Cluster")
        axes[1].set_xticklabels(dictionary.keys(), rotation = 90)
        plt.margins(x=0.01)
        fig.tight_layout()
    else: 
        fig, ax = plt.subplots(figsize = (width, height))
        title = f"{name_E1} ({len(np.unique(adata_E1.obs[cluster_header_E1]))} clusters, {adata_E1.shape[0]} cells)"
        dictionary = dict(adata_E1.obs[cluster_header_E1].value_counts())
        if proportion: 
            total = sum(dictionary.values())
            for key in dictionary: dictionary[key] =round(dictionary[key]/total, 3)
        p = ax.bar(dictionary.keys(), dictionary.values(), color = color_E1)
        ax.set_ylim(0, list(dictionary.values())[0] + y_offset)
        ax.bar_label(p, fontsize = 8, rotation = 30)
        ax.set_title(title)
        if proportion: ax.set_ylabel("Proportion")
        else: ax.set_ylabel("Size")
        ax.set_xlabel("Cluster")
        ax.set_xticklabels(dictionary.keys(), rotation = 90)
        plt.margins(x=0.01)

    if save == True: 
        plt.savefig(f"clusterSize.png")
    elif isinstance(save, str): 
        plt.savefig(save)
    return
    