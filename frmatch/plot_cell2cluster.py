
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

def plot_cell2cluster(results, reorder = True, order_query = [], order_ref = [], axis = 0, p_adj_method = "fdr_by", sig_level = 0.1, prefix = ["query_", "ref_"], width = 10, height = 10, title = None, save_intermediate = False, save = False): 

    folder = "tutorial_hlca_cellref"
    # Pivotting results to ref_cluster as index and query_cluster as columns
    results = results.sort_values(["query_cluster", "index", "p_value"], ascending = [True, True, False])
    results = results[~results.duplicated(["index", "query_cluster", "ref_cluster"], keep = "first")]
    results = results.pivot(index = "ref_cluster", columns = ["query_cluster", "index"], values = "p_value").replace(np.nan, 0)
    if save_intermediate: results.to_csv(f"results_{p_adj_method}_{sig_level}.csv")

    results_2 = FRmatch.padj_FRmatch(results, p_adj_method = p_adj_method)
    if save_intermediate: results_2.to_csv(f"results_2_{p_adj_method}_{sig_level}.csv")
    
    # Setting values to zero unless passes sig_level threshold
    results_3 = results_2.applymap(lambda x: x if x > sig_level else 0)
    if save_intermediate: results_3.to_csv(f"results_3_{p_adj_method}_{sig_level}.csv")

    # Adding the unassigned row
    samples = results_3.columns[results_3.sum() == 0]
    temp = pd.DataFrame(dict(zip(samples, [1] * len(samples))), index = ["unassigned"])
    if temp.shape[1] != 0: 
        temp.columns.names = ["query_cluster", "index"]
        temp.index.name = "ref_cluster"
        results_4 = pd.concat([results_3, temp])
        if save_intermediate: results_4.to_csv(f"results_4_{p_adj_method}_{sig_level}.csv")
    else: 
        results_4 = results_3.copy()
    
    # idxmax: returns ref_cluster with largest value per sample
    results_5 = pd.DataFrame(results_4.idxmax()).reset_index().rename(columns = {0: "ref_cluster"})
    results_5 = pd.DataFrame(results_5.value_counts(["query_cluster", "ref_cluster"])).reset_index()
    results_5 = results_5.pivot(index='ref_cluster', columns='query_cluster')
    if save_intermediate: results_5.to_csv(f"results_5_{p_adj_method}_{sig_level}.csv")
    
    results_6 = results_5.replace(np.nan, 0)
    results_6.columns = [val[1] for val in results_6.columns]
    results_6.columns.name = "query_cluster"
    if reorder: 
        # ensuring unassigned row is at the bottom
        unassigned = results_6.loc[results_6.index == 'unassigned'].copy()
        results_6 = results_6.loc[results_6.index != 'unassigned'].copy()
        results_6 = FRmatch.reorder_FRmatch(results_6, axis = axis)
        results_6 = pd.concat([results_6, unassigned])
    if save_intermediate: results_6.to_csv(f"results_6_{p_adj_method}_{sig_level}.csv")
    
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
    if save_intermediate: results_7.to_csv(f"results_7_{p_adj_method}_{sig_level}.csv")
    
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

    if save == True: 
        plt.savefig(f"frmatch_results_cell2cluster.png")
    elif isinstance(save, str): 
        plt.savefig(save)
    
    return

# def plot_cell2cluster_old(results, reorder = True, order_query = [], order_ref = [], axis = 0, p_adj_method = "fdr_by", sig_level = 0.1, prefix = ["query_", "ref_"], width = 6, height = 6, title = None, save = False): 

#     # Sorting my query_cluster per sample's (index) p_value
#     results = results.sort_values(["query_cluster", "index", "p_value"], ascending = [True, True, False])
#     # Checking for duplicated samples (index)
#     results = results[~results.duplicated(["index", "query_cluster", "ref_cluster"], keep = "first")]    # Grouping all samples (index) by query_clusters
#     results = results.pivot(index = "ref_cluster", columns = ["query_cluster", "index"], values = "p_value").replace(np.nan, 0)
#     results = FRmatch.padj_FRmatch(results, p_adj_method = p_adj_method)
#     zeros = []
#     for col in results.columns: 
#         maximum = max(results[col])
#         if maximum < sig_level: 
#             maximum = -1
#         values = results[col] == maximum
#         results[col] = values
#         if sum(results[col]) == 0: 
#             zeros.append(col)
#     temp = pd.DataFrame(dict(zip(results.columns, [1 if val in zeros else 0 for val in results.columns])), index = ["unassigned"])
#     temp.columns.names = ["query_cluster", "index"]
#     results = pd.concat([results, temp])
#     results.index.name = "ref_cluster"
#     results.to_csv("results.csv")
    
#     results_2 = results.groupby(level = ["query_cluster"], axis = 1).sum()
#     if reorder: 
#         unassigned = results_2.loc[results_2.index == 'unassigned'].copy()
#         results_2 = results_2.loc[results_2.index != 'unassigned'].copy()
#         results_2 = FRmatch.reorder_FRmatch(results_2, axis = axis)
#         results_2 = pd.concat([results_2, unassigned])
#     results_2.to_csv("results_2.csv")
#     results_3 = pd.DataFrame(columns = ["ref_cluster", "query_cluster", "value"])
#     results_3 = pd.DataFrame(columns = ["match", "query_cluster", "Prop"])
#     for index, row in results_2.iterrows(): 
#         for col in results_2.columns: 
#             if sum(results_2[col]) == 0: value = 0
#             else: value = row[col] / sum(results_2[col])
#             results_3 = pd.concat([results_3, pd.DataFrame(dict(zip(results_3.columns, [index, col, value])), index = [0])])
    
#     if not reorder and len(order_query) != 0: 
#         # query cluster order
#         results_3["query_cluster"] = results_3["query_cluster"].astype("category")
#         results_3["query_cluster"] = results_3["query_cluster"].cat.set_categories(order_query)
#     if not reorder and len(order_ref) != 0: 
#         # reference cluster order
#         results_3["match"] = results_3["match"].astype("category")
#         results_3["match"] = results_3["match"].cat.set_categories(order_ref)
#         results_3 = results_3.sort_values(["match", "query_cluster"])
#     results_3["match"] = [prefix[0] + val for val in results_3["match"]]
#     results_3["query_cluster"] = [prefix[1] + val for val in results_3["query_cluster"]]
#     results_3.to_csv("results_3.csv")
#     if not title: title = "FR-Match cell-to-cluster"
#     fig, (ax1) = plt.subplots(1, 1, figsize=(width, height))
#     sns.set_theme(style="whitegrid", rc={"grid.color": "lightgray", "grid.linewidth": 0.5})
#     # sns.set_theme(style="whitegrid", rc = {"axes.edgecolor": "darkgray", "xtick.bottom": True, "ytick.left": True})
#     ax = sns.scatterplot(results_3, x = "query_cluster", y = "match", hue = "Prop", size = "Prop", 
#                         hue_norm = (0, 1), size_norm = (0, 1), legend = "brief", 
#                         sizes = (1, 200), alpha = 0.75, linewidth = 0.7, edgecolor = "gray", 
#                         palette = sns.color_palette("viridis", as_cmap=True), 
#                        )
#     a = plt.title(title)
#     a = plt.xticks(rotation = 90)
#     a = plt.legend(title = "Prop", bbox_to_anchor = (1.0, 0.7))

#     if save == True: 
#         plt.savefig(f"frmatch_results_cell2cluster.png")
#     elif isinstance(save, str): 
#         plt.savefig(save)
    
#     return 
