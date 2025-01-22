
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib

def plot_clusterSize(adata_E1, adata_E2, cluster_header, name_E1 = "E1", name_E2 = "E2", 
                     color_E1 = "#F0E442", color_E2 = "#56B4E9", width = 10, height = 10): 
    
    fig, axes = plt.subplots(2, figsize = (width, height))
    
    title = f"{name_E1} ({len(np.unique(adata_E1.obs[cluster_header]))} clusters, {adata_E1.shape[0]} cells)"
    dictionary = dict(adata_E1.obs[cluster_header].value_counts())
    p = axes[0].bar(dictionary.keys(), dictionary.values(), color = color_E1)
    axes[0].bar_label(p, fontsize = 8)
    axes[0].set_title(title)
    axes[0].set_ylabel("Size")
    axes[0].set_xlabel("Cluster")
    axes[0].set_xticklabels(dictionary.keys(), rotation = 90)
    
    title = f"{name_E2} ({len(np.unique(adata_E2.obs[cluster_header]))} clusters, {adata_E2.shape[0]} cells)"
    dictionary = dict(adata_E2.obs[cluster_header].value_counts())
    p = axes[1].bar(dictionary.keys(), dictionary.values(), color = color_E2)
    axes[1].bar_label(p, fontsize = 8)
    axes[1].set_title(title)
    axes[1].set_ylabel("Size")
    axes[1].set_xlabel("Cluster")
    axes[1].set_xticklabels(dictionary.keys(), rotation = 90)
    
    fig.tight_layout()
    
    return
