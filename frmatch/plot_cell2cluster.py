
import os
import sys
sys.path.insert(0, os.path.abspath("./frmatch"))
import frmatch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from random import sample
import matplotlib.patches as mpatches

def plot_cell2cluster(results, reorder = True, order_query = [], order_ref = [], axis = 0, p_adj_method = "fdr_by", sig_level = 0.1, prefix = ["query_", "ref_"], width = 10, height = 10, title = None, save = False): 
    results_5 = pd.DataFrame(results.value_counts(["query_cluster", "ref_cluster"])).reset_index().rename(columns = {0: "count"})
    results_5 = results_5.pivot(index='ref_cluster', columns='query_cluster').replace(np.nan, 0)
    results_5.columns = [val[1] for val in list(results_5.columns)]
    results_5.columns.name = "query_cluster"
    
    if reorder: 
        # ensuring unassigned row is at the bottom
        unassigned = results_5.loc[results_5.index == 'unassigned'].copy()
        results_6 = results_5.loc[results_5.index != 'unassigned'].copy()
        results_6 = frmatch.reorder_FRmatch(results_6, axis = axis)
        results_6 = pd.concat([results_6, unassigned])
    else: results_6 = results_5.copy()
    
    results_7 = pd.DataFrame(columns = ["match", "query_cluster", "Prop"]) # ref_cluster, query_cluster, proportion
    for index, row in results_6.iterrows(): 
        for col in results_6.columns: 
            if sum(results_6[col]) == 0: value = 0
            else: value = row[col] / sum(results_6[col])
            results_7 = pd.concat([results_7, pd.DataFrame(dict(zip(results_7.columns, [index, col, value])), index = [0])])
    if not reorder and len(order_query) != 0: 
        # query cluster order
        results_7["query_cluster"] = results_7["query_cluster"].astype("category")
        results_7["query_cluster"] = results_7["query_cluster"].cat.set_categories(order_query)
    if not reorder and len(order_ref) != 0: 
        # reference cluster order
        results_7["match"] = results_7["match"].astype("category")
        results_7["match"] = results_7["match"].cat.set_categories(order_ref)
        results_7 = results_7.sort_values(["match", "query_cluster"])
    results_7["match"] = [prefix[0] + val for val in results_7["match"]]
    results_7["query_cluster"] = [prefix[1] + val for val in results_7["query_cluster"]]
    
    if not title: title = "FR-Match cell-to-cluster"
    fig, (ax1) = plt.subplots(1, 1, figsize=(width, height))
    sns.set_theme(style="whitegrid", rc={"grid.color": "lightgray", "grid.linewidth": 0.5})
    ax = sns.scatterplot(results_7, x = "query_cluster", y = "match", hue = "Prop", size = "Prop", 
                        hue_norm = (0, 1), size_norm = (0, 1), legend = "brief", 
                        sizes = (1, 200), alpha = 0.75, linewidth = 0.7, edgecolor = "gray", 
                        palette = sns.color_palette("viridis", as_cmap=True), 
                       )
    a = plt.title(title)
    a = plt.xticks(rotation = 90)
    a = plt.legend(title = "Prop", bbox_to_anchor = (1.0, 0.7))
    sns.set_theme(style="whitegrid", rc={"grid.color": "lightgray", "grid.linewidth": 0.5})

    if save == True: 
        plt.savefig(f"frmatch_results_cell2cluster.png")
    elif isinstance(save, str): 
        plt.savefig(save)
    
    return
