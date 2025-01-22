
from .utils import get_markers
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from random import sample
# from pheatmap import pheatmap
import matplotlib.patches as mpatches

def plot_cluster_by_markers(query, cluster_header, cluster, markers = None, marker_legend_loc = None, 
                            nsamp = 30, colorbar_loc = (0, 0.88), 
                            name_adata = "E1", name_markers = "E2", use_common_markergenes = True,
                            scale_colorbar = False, cellheight = 10, cellwidth = 5, title = None, filename = None): 
    """\
    
    """
    ## check if the cluster is in query
    if cluster not in list(np.unique(query.obs[cluster_header])): 
        return f"{cluster} not found in query"
    
    ## marker genes IN ORDER
    if markers: 
        marker_genes = markers[:]
    else: 
        marker_genes = get_markers(query, cluster_header, "NSForest_markers")
    if not title: 
        title = f"{name_markers} markers on {name_adata}\n{cluster}"
#     print(title)
    ## cells of query cluster
    col_query = query.obs[cluster_header] == cluster
    
    if use_common_markergenes: 
        marker_genes_common = []
        for marker in marker_genes: 
            if marker in list(query.var_names): 
                marker_genes_common.append(marker)
        marker_genes = marker_genes_common[:]
    query_df = query[query.obs_names[col_query], marker_genes].to_df()
    
    ## randomly select nsamp number of cells
    if query_df.shape[0] > nsamp: 
        mm = sample(range(query_df.shape[0]), nsamp)
        mm = [list(query_df.index)[val] for val in mm]
        query_df = query_df.loc[mm,:]
    
    query_df = query_df.T
    if query_df.shape[1] < 30: print(query_df.shape)
    
    ## if to scale colorbar to [0,1] for normalized expr

    # Plotting 
    fig, (ax1) = plt.subplots(1, 1, figsize = (6, 12))
    ax = sns.heatmap(query_df, cmap = "inferno", vmin = 0, cbar = False, yticklabels = 1, square = True, ax = ax1)
    if len(title) < 26: a = plt.title(title, size = 12)
    else: a = plt.title(title, size = 10)
    a = plt.ylabel("")
    ax.yaxis.tick_right()
    a = plt.yticks(rotation = 0, size = 6)
    a = plt.xticks([], [])
    pcm = ax1.pcolormesh(query_df, cmap = "inferno", vmin = 0)
    a = fig.colorbar(pcm, ax=ax1, shrink=0.2, location = "right", pad = 0.15, label = "expression (log2norm)", 
                     anchor = colorbar_loc) # anchor = (0, 0.5)

    if marker_legend_loc: 
        cluster_markers = list(query.uns["nsforest_results"][cluster_header][cluster]["NSForest_markers"])
        row_colors = ["#00BFC4" if marker in cluster_markers else "lightgray" for marker in marker_genes]
        for i, color in enumerate(row_colors):
            ax.add_patch(plt.Rectangle(xy=(-0.07, i), width=0.05, height=1, color=color, lw=0, transform=ax.get_yaxis_transform(), clip_on=False))
        handles = [mpatches.Patch(color='lightgray', label='0'), 
                   mpatches.Patch(color='#00BFC4', label='1')]
        a = plt.legend(title = "Marker", handles = handles, bbox_to_anchor = marker_legend_loc)

    if filename: 
        plt.savefig(filename, bbox_inches = "tight")
        plt.close()
    return
