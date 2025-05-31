
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

def plot_cell2cluster(results, reorder = True, order_query = [], order_ref = [], axis = 0, sig_level = 0.1, prefix = ["query", "ref"]): 
    # Sorting my query_cluster per sample's (index) p_value
    results = results.sort_values(["query_cluster", "index", "p_value"], ascending = [True, True, False])
    # Checking for duplicated samples (index)
    results = results[~results.duplicated(["index", "query_cluster", "ref_cluster"], keep = "first")]    # Grouping all samples (index) by query_clusters
    results = results.pivot(index = "ref_cluster", columns = ["query_cluster", "index"], values = "p_value").replace(np.nan, 0)
    results = FRmatch.padj_FRmatch(results, p_adj_method = "fdr_bh")
    zeros = []
    for col in results.columns: 
        maximum = max(results[col])
        if maximum < sig_level: 
            maximum = -1
        values = results[col] == maximum
        results[col] = values
        if sum(results[col]) == 0: 
            zeros.append(col)
    temp = pd.DataFrame(dict(zip(results.columns, [1 if val in zeros else 0 for val in results.columns])), index = ["unassigned"])
    temp.columns.names = ["query_cluster", "index"]
    results = pd.concat([results, temp])
    results.index.name = "ref_cluster"
    
    results_2 = results.groupby(level = ["query_cluster"], axis = 1).sum()
    if reorder: 
        unassigned = results_2.loc[results_2.index == 'unassigned'].copy()
        results_2 = results_2.loc[results_2.index != 'unassigned'].copy()
        results_2 = FRmatch.reorder_FRmatch(results_2, axis = axis)
        results_2 = pd.concat([results_2, unassigned])
    
    results_3 = pd.DataFrame(columns = ["ref_cluster", "query_cluster", "value"])
    results_3 = pd.DataFrame(columns = ["match", "query_cluster", "Prop"])
    for index, row in results_2.iterrows(): 
        for col in results_2.columns: 
            if sum(results_2[col]) == 0: value = 0
            else: value = row[col] / sum(results_2[col])
            results_3 = pd.concat([results_3, pd.DataFrame(dict(zip(results_3.columns, [index, col, value])), index = [0])])
    
    if not reorder and len(order_query) != 0: 
        # query cluster order
        results_3["query_cluster"] = results_3["query_cluster"].astype("category")
        results_3["query_cluster"] = results_3["query_cluster"].cat.set_categories(order_query)
    if not reorder and len(order_ref) != 0: 
        # reference cluster order
        results_3["match"] = results_3["match"].astype("category")
        results_3["match"] = results_3["match"].cat.set_categories(order_ref)
        results_3 = results_3.sort_values(["match", "query_cluster"])
    results_3["match"] = [prefix[0] + val for val in results_3["match"]]
    results_3["query_cluster"] = [prefix[1] + val for val in results_3["query_cluster"]]

    fig, (ax1) = plt.subplots(1, 1, figsize=(6, 6))
    sns.set_theme(style="whitegrid", rc={"grid.color": "lightgray", "grid.linewidth": 0.5})
    # sns.set_theme(style="whitegrid", rc = {"axes.edgecolor": "darkgray", "xtick.bottom": True, "ytick.left": True})
    ax = sns.scatterplot(results_3, x = "query_cluster", y = "match", hue = "Prop", size = "Prop", 
                        hue_norm = (0, 1), size_norm = (0, 1), legend = "brief", 
                        sizes = (1, 200), alpha = 0.75, linewidth = 0.7, edgecolor = "gray", 
                        palette = sns.color_palette("viridis", as_cmap=True), 
                       )
    a = plt.xticks(rotation = 90)
    a = plt.legend(title = "Prop", bbox_to_anchor = (1.0, 0.7))
    
    return 
