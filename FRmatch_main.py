
import os
import sys
sys.path.insert(0, os.path.abspath("./FRmatch"))
import FRmatch
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm

def FRmatch_main(query, ref, cluster_header, marker_col = "markers", use_cosine = True, imputation = False, 
            filter_size = 5, filter_fscore = None, filter_nomarker = False, 
            add_pseudo_marker = False, pseudo_expr = 1, 
            subsamp_size = 20, subsamp_iter = 1000, subsamp_seed = 916, 
            numCores = None, prefix = ["query", "ref"], 
            verbose = 0, return_all = False, save_as = False): 
    """\
    
    """
    if save_as: 
        if not isinstance(save_as, str): 
            filename = f"frmatch_results_{prefix[0]}to{prefix[1]}.csv"
        elif ".csv" not in str(save_as): 
            filename = filename + ".csv"
        else: 
            filename = save_as
        
    ## check data object

    ## find common marker genes
    marker_genes = FRmatch.get_markers(ref, cluster_header, marker_col)
    
    marker_genes_common = []
    for marker in marker_genes: 
        if marker in list(query.var_names): 
            marker_genes_common.append(marker)
    if verbose > 0: 
        print(f"Feature selection: {len(marker_genes_common)} out of {len(marker_genes)} reference marker genes are presented in the query experiment.")
    if len(marker_genes) != len(marker_genes_common) and verbose > 0: 
        print(f"{sorted(list(set(marker_genes) - set(marker_genes_common)))} not found.")
    
    ## filtering ref clusters without marker genes available in query
    
    ## filtering small or low fscore clusters
    if verbose > 0: 
        print(f"Filtering small clusters: query and reference clusters with less than {filter_size} cells are not considered.")
    if filter_fscore != None: 
        return
    
    query_shape = query.shape
    ref_shape = ref.shape
    query = FRmatch.filter_cluster(query, cluster_header, filter_size) #only filter on size, not fscore
    ref = FRmatch.filter_cluster(ref, cluster_header, filter_size)
    if verbose > 0: 
        print(f"filtered query from {query_shape} -> {query.shape}")
        print(f"filtered ref from {ref_shape} -> {ref.shape}")
    
    ## get expr data and dimension reduction by selecting common marker genes
    query_X_filt = query[:, marker_genes_common].to_df()
    query_X_filt["cluster"] = query.obs[cluster_header]
    ref_X_filt = ref[:, marker_genes_common].to_df()
    ref_X_filt["cluster"] = ref.obs[cluster_header]
    
    ## if to add pseudo marker to give signals in query clusters with almost no expression
    
    ## extract cluster info from sce.objects
    query_clusters = list(query.obs[cluster_header])
    ref_clusters = list(ref.obs[cluster_header])

    ## split data by cluster membership and store in a list
    
    ## FR comparison between two experiments
    ## prepare data for each combination of query cluster and ref cluster pair
    
    results = pd.DataFrame()
    for query_cluster in tqdm(np.unique(query_clusters), desc = "FRmatch"): 
        if verbose > 0: print("QUERY_CLUSTER: ", query_cluster)

        for ref_cluster in np.unique(ref_clusters): 
            if verbose > 0: print("ref_cluster: ", ref_cluster)
            
            # subsetting query and ref to query_cluster and ref_cluster
            query_df = query_X_filt[query_X_filt["cluster"] == query_cluster]
            del query_df["cluster"]
            ref_df = ref_X_filt[ref_X_filt["cluster"] == ref_cluster]
            del ref_df["cluster"]
            df = FRmatch.FRtest_subsamp(query_df, ref_df, subsamp_size = subsamp_size, subsamp_iter = subsamp_iter, subsamp_seed = subsamp_seed, return_all = return_all)
            df["query_cluster"] = query_cluster
            df["ref_cluster"] = ref_cluster
            results = pd.concat([results, df])

        if save_as: 
            results_temp = results.pivot(index = "ref_cluster", columns = "query_cluster", values = "p_value")
            results_temp.to_csv(filename, index = True)
        
    results_pval = results.pivot(index = "ref_cluster", columns = "query_cluster", values = "p_value")
    if save_as: 
        print(f"Saving p_values as... {filename}")
        results_pval.to_csv(filename, index = True)

    return results_pval
