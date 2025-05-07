
import os
import sys
sys.path.insert(0, os.path.abspath("./FRmatch"))
import FRmatch

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import matplotlib.patches as mpatches

def plot_FRmatch(pmat, type_ = "matches", p_adj_method="fdr_by", sig_level = 0.05, marker_legend_loc = (2.4, 1), 
                 reorder = True, ignore_unassigned = False, return_value = False, colors = ['#4575B4', '#FEE090'], 
                 figsize_width = 6, figsize_height = 12, main = None): 
    
    ## calculate adjusted p-values and determine matches
    pmat_adj = FRmatch.padj_FRmatch(pmat, p_adj_method = p_adj_method)
    pmat_cutoff = FRmatch.cutoff_FRmatch(pmat, p_adj_method = p_adj_method, sig_level = sig_level)
    
    ## reorder
    if reorder: 
        pmat_cutoff = FRmatch.reorder_FRmatch(pmat_cutoff)
        pmat_adj = pmat_adj[list(pmat_cutoff.columns)]
    
    ## plot
    if type_ == "matches": 
        if not main: main = "FR-Match cluster-to-cluster"
        if ignore_unassigned: 
            pmat_cutoff = pmat_cutoff.drop("unassigned")
        fig, (ax1) = plt.subplots(1, 1, figsize=(figsize_width, figsize_height))
        ax = sns.heatmap(pmat_cutoff, cmap = colors, cbar = False, yticklabels = 1, square = True, ax = ax1, 
                         linewidths = 0.5, linecolor = "gray") # cmap = "RdYlBu", 
        a = plt.title(f"{main}")
        ax.yaxis.tick_right()
        a = plt.yticks(rotation = 0, size = 6)
        a = plt.xticks(rotation = 270, size = 6)
        if marker_legend_loc: 
            handles = [mpatches.Patch(color=colors[1], label='Match'), 
                       mpatches.Patch(color=colors[0], label='No match')] 
            a = plt.legend(title = "", handles = handles, bbox_to_anchor = marker_legend_loc) # (1.53, 1) (2, 1)
            
    elif type_ == "padj": 
        values = pd.DataFrame(pmat_adj.unstack()).reset_index()
        fig, (ax1) = plt.subplots(1, 1, figsize=(figsize_width, figsize_height))
        a = plt.scatter(values["level_0"], values[0], color = "black")
        a = plt.title(f"{main} {p_adj_method} {sig_level}")
        a = plt.ylabel("Adjusted p-value")
        a = plt.xlabel("Query cluster")
        a = plt.xticks(rotation = 270)
        a = plt.grid(alpha = 0.25)
        a = plt.axhline(sig_level, color = "red", linestyle = '--', alpha = 0.5)
    
    return
