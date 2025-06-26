
import os
import sys
sys.path.insert(0, os.path.abspath("./FRmatch"))
import FRmatch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib
import seaborn as sns

def plot_FRmatch(pmat, type_ = "matches", p_adj_method="fdr_by", sig_level = 0.05, marker_legend_loc = (2.4, 1), 
                 reorder = True, ignore_unassigned = False, width = 6, height = 12, title = None, save = False): 
    """\
    Plotting one-directional FRmatch matching results after running padj_FRmatch and cutoff_FRmatch.

    Parameters
    ----------
        pmat: pd.DataFrame
            P-value matrix.
        type_: ["matches", "padj"] (default: "matches")
            Type of plot. If "matches", plots heatmap of cutoff_FRmatch. If "padj", plots adjusted p-values with `sig_level` as theshold.
        p_adj_method: str (default: "fdr_by")
            P-value adjustment method. 
        sig_level: float (default: 0.05)
            P-value cutoff threshold.
        marker_legend_loc: tuple (default: (2.4, 1))
            Location of legend.
        reorder: bool (default: True)
            Whether to run reorder_FRmatch into a diagonal.
        ignore_unassigned: bool (default: False)
            Whether remove the unassigned row from "matches" heatmap.
        width: int (default: 6)
            Width of plot.
        height: int (default: 12)
            Height of plot.
        title: str (default: None)
            Plot title.
    """
    ## calculate adjusted p-values and determine matches
    pmat_adj = FRmatch.padj_FRmatch(pmat, p_adj_method = p_adj_method)
    pmat_cutoff = FRmatch.cutoff_FRmatch(pmat, p_adj_method = p_adj_method, sig_level = sig_level)
    
    ## reorder
    if reorder: 
        pmat_cutoff = FRmatch.reorder_FRmatch(pmat_cutoff)
        pmat_adj = pmat_adj[list(pmat_cutoff.columns)]
    
    ## plot
    if type_ == "matches": 
        if not title: title = "FR-Match cluster-to-cluster"
        if ignore_unassigned: 
            pmat_cutoff = pmat_cutoff.drop("unassigned")
        fig, (ax1) = plt.subplots(1, 1, figsize=(width, height))
        ax = sns.heatmap(pmat_cutoff, cmap = ['#4575B4', '#FEE090'], cbar = False, yticklabels = 1, 
                         square = True, ax = ax1, linewidths = 0.5, linecolor = "gray") 
        a = plt.title(f"{title}")
        ax.yaxis.tick_right()
        a = plt.yticks(rotation = 0) # size = 6
        a = plt.xticks(rotation = 270) # size = 6
        if marker_legend_loc: 
            handles = [mpatches.Patch(color='#FEE090', label='Match'), 
                       mpatches.Patch(color='#4575B4', label='No match')]
            a = plt.legend(title = "", handles = handles, bbox_to_anchor = marker_legend_loc)
            
    elif type_ == "padj": 
        if not title: title = "FR-Match cluster-to-cluster adjusted p-values"
        values = pd.DataFrame(pmat_adj.unstack()).reset_index()
        fig, (ax1) = plt.subplots(1, 1, figsize=(width, height))
        a = plt.scatter(values["level_0"], values[0], color = "black")
        a = plt.title(f"{title}")
        a = plt.ylabel("Adjusted p-value")
        a = plt.xlabel("Query cluster")
        a = plt.xticks(rotation = 270)
        a = plt.grid(alpha = 0.25)
        a = plt.axhline(sig_level, color = "red", linestyle = '--', alpha = 0.5)

    if save == True: 
        plt.savefig(f"frmatch_results_onedirectional_{type_}.png")
    elif isinstance(save, str): 
        plt.savefig(save)
    
    return
