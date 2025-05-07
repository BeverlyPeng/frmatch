
import os
import sys
sys.path.insert(0, os.path.abspath("./FRmatch"))
import FRmatch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from random import sample
import matplotlib.patches as mpatches

def plot_cell2cluster(e1_e2, e2_e1, prefix = ["query", "ref"], axis = 0,
                      p_adj_method="fdr_by", sig_level = 0.05, marker_legend_loc = (2.6, 1), 
                      reorder = True, two_way_only = False, return_value = False,
                      cellwidth = 10, cellheight = 10, title = None, filename = ""): 
    
    ## get binary matrices for plotting
    pmat_cutoff_e1_e2 = FRmatch.cutoff_FRmatch(e1_e2, p_adj_method = p_adj_method, sig_level = sig_level)
    pmat_cutoff_e2_e1 = FRmatch.cutoff_FRmatch(e2_e1, p_adj_method = p_adj_method, sig_level = sig_level)
    
    ## combine two matrices to one two-way matrix
    mat1 = pmat_cutoff_e1_e2.drop("unassigned")
    mat2 = pmat_cutoff_e2_e1.drop("unassigned").T
    mat_bi = mat1 + mat2
    
    if two_way_only: 
        mat_bi
    
    ## unassigned row
    mat_bi = pd.concat([mat_bi, pd.DataFrame(dict(2*(mat_bi.sum() == 0)), index = ["unassigned"])])
        
#     ## rename colnames and rownames
#     mat_bi.index = [prefix[1] + name_e2 + val for val in pmat_cutoff_e1_e2.index]
#     mat_bi.columns = [prefix[0] + name_e1 + val for val in pmat_cutoff_e1_e2.columns]
    
    ## plot
    if not title: title = "FR-Match cluster-to-cluster"
    if reorder: 
        mat_bi = FRmatch.reorder_FRmatch(mat_bi, axis = axis)
    if two_way_only: 
        mat_bi
    else: 
        fig, (ax1) = plt.subplots(1, 1, figsize=(12, 12))
        ax = sns.heatmap(mat_bi, cmap = ["#4575B4", "#FEE090", "#D73027"], cbar = False, 
                         yticklabels = 1, square = True, ax = ax1, 
                         linewidths = 0.5, linecolor = "gray") # cmap = "RdYlBu", 
        a = plt.title(f"{title}")
        a = plt.xlabel(f"{prefix[0]}clusters")
        a = plt.ylabel(f"{prefix[1]}clusters")
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        a = plt.yticks(rotation = 0)
        a = plt.xticks(rotation = 270)
        if marker_legend_loc: 
            handles = [mpatches.Patch(color='#D73027', label='Two-way match'), 
                       mpatches.Patch(color='#FEE090', label='One-way match'), 
                       mpatches.Patch(color='#4575B4', label='No match')] 
            a = plt.legend(title = "", handles = handles, bbox_to_anchor = marker_legend_loc) # (1.53, 1)
    
    return 
