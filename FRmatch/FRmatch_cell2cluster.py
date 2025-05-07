
from IPython.core.debugger import set_trace
import os
import sys
sys.path.insert(0, os.path.abspath("./FRmatch"))
import FRmatch
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm
import pickle

def FRmatch_cell2cluster(query, ref, cluster_header_query, cluster_header_ref, use_cosine = True, #imputation=FALSE,
                         feature_selection = "reference_markers", #feature selection
                         filter_size = 5, filter_fscore = None, filter_nomarker = False, #filtering clusters
                         add_pseudo_marker = False, pseudo_expr = 1, #adding pseudo marker
                         subsamp_size = 10, subsamp_iter = 2000, subsamp_seed = 916, #subsampling
                         subsamp_iter_custom = False, subsamp_iter_custom_k = 5, #customization
                         numCores = None, prefix = ["query", "ref"], verbose = 0, save_as = False): 
    
    settings = {"query": prefix[0], "ref": prefix[1], "cluster_header_query": cluster_header_query, "cluster_header_ref": cluster_header_ref, "feature_selection": feature_selection, "use_cosine": use_cosine, "filter_size": filter_size, "subsamp_size": subsamp_size, "subsamp_iter": subsamp_iter, "subsamp_seed": subsamp_seed, "subsamp_iter_custom": subsamp_iter_custom, "subsamp_iter_custom_k": subsamp_iter_custom_k, "save_as": save_as}
    
    if save_as: 
        if not isinstance(save_as, str): 
            filename = f"frmatch_results_cell2cluster_{prefix[0]}to{prefix[1]}.pkl"
        elif ".pkl" not in str(save_as): 
            filename = save_as + ".pkl"
        else: 
            filename = save_as

    ## check data object
    
    ## feature selection
    if feature_selection == "query_genes": 
        marker_genes = list(query.var_names)
        marker_genes_common = list(set(marker_genes) & set(ref.var_names))
        if verbose > 0: 
            print(f"Feature selection: marker genes of the reference experiment.")
            print(f"\t{len(marker_genes_common)} out of {len(marker_genes)} query genes are common in the reference experiment.")
        if len(marker_genes) != len(marker_genes_common) and verbose > 0: 
            print(f"{sorted(list(set(marker_genes) - set(marker_genes_common)))} not found.")
    elif feature_selection == "reference_markers": 
#         marker_genes = ref.uns["markers"]
        marker_genes = list(ref.var_names)
        marker_genes_common = list(set(marker_genes) & set(query.var_names))
        if verbose > 0: 
            print(f"Feature selection: gene space of the query experiment.")
            print(f"\t{len(marker_genes_common)} out of {len(marker_genes)} reference marker genes are presented in the query experiment.")
        if len(marker_genes) != len(marker_genes_common) and verbose > 0: 
            print(f"{sorted(list(set(marker_genes) - set(marker_genes_common)))} not found.")
    elif feature_selection == "all": 
        marker_genes_common = list(set(query.var_names) & set(ref.var_names))
        if verbose > 0: 
            print(f"Feature selection: intersection of genes available from query and ref.")
    ## filtering ref clusters without marker genes available in query
    
    ## filtering small or low fscore clusters
#     if verbose > 0: 
#         print(f"Filtering small clusters: query and reference clusters with less than {filter_size} cells are not considered.")
#     if filter_fscore != None: 
#         return
        
    ## reduced data
    query_reduced = query[:, marker_genes_common]
    ref_reduced = ref[:, marker_genes_common]
    print(query_reduced.shape)
    print(ref_reduced.shape)
    
    # Decide which step later
    query_X_filt = query[:, marker_genes_common].to_df()
    query_X_filt["cluster"] = query.obs[cluster_header_query]
    ref_X_filt = ref[:, marker_genes_common].to_df()
    ref_X_filt["cluster"] = ref.obs[cluster_header_ref]
    
    ## add pseudo marker
    
    ## data prep
    ## extract cluster info from sce.objects
    membership_query = query_reduced.obs[cluster_header_query]
    membership_ref = ref_reduced.obs[cluster_header_ref]
    # order_query = list(query_reduced.uns[f"dendrogram_{cluster_header_query}"]["categories_ordered"])
    # order_ref = list(ref_reduced.uns[f"dendrogram_{cluster_header_ref}"]["categories_ordered"])

    ## split data by cluster membership and store in a list
    
    ncluster_query = len(np.unique(membership_query))
    ncluster_ref = len(np.unique(membership_ref))
    if verbose > 0: 
        print(f"Comparing {ncluster_query} query clusters with {ncluster_ref} reference clusters.")
    
    ## FR comparison between two experiments
    results = pd.DataFrame()
    for query_cluster in tqdm(np.unique(membership_query), desc = "FRmatch"): 
        if verbose > 0: print("QUERY_CLUSTER:", query_cluster)
        warning_message = True
        for ref_cluster in np.unique(membership_ref): 
            if verbose > 0: print("ref_cluster:", ref_cluster)
            
            # subsetting query and ref to query_cluster and ref_cluster
            query_df = query_X_filt[query_X_filt["cluster"] == query_cluster]
            del query_df["cluster"]
            ref_df = ref_X_filt[ref_X_filt["cluster"] == ref_cluster]
            del ref_df["cluster"]
#             print(query_df.shape, ref_df.shape)
            
            if subsamp_iter_custom: 
                subsamp_iter = max(subsamp_iter, subsamp_iter_custom_k * query_df.shape[0])
                
            df = FRmatch.FRtest_cell2cluster(query_df, ref_df, subsamp_size = subsamp_size, subsamp_iter = subsamp_iter, subsamp_seed = subsamp_seed)
            
            if df.shape[0] != query_df.shape[0] and warning_message: 
                print(f"warning: Not all cells are randomly sampled. Consider increasing the number of iterations specified in the 'subsamp_iter' argument. query_cluster: {query_cluster} ref_cluster: {ref_cluster}")
                warning_message = False
            
            df["query_cluster"] = query_cluster
            df["ref_cluster"] = ref_cluster
            results = pd.concat([results, df])
            if save_as: 
#                 print(f"Saving p_values as... {filename}")
#                 results.to_csv(filename, index = True)
                with open(filename, 'wb') as f:
                    pickle.dump({'settings': settings, "results": results}, f)
            
#             if(subsamp.iter.custom) subsamp.iter = max(subsamp.iter, subsamp.iter.custom.k*ncol(samp1))
#               ## start iterations by subsampling
#               set.seed(subsamp.seed)
#               FRtest_cell2cluster(samp1, samp2, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, ...)
    if save_as: 
#         print(f"Saving p_values as... {filename}")
#         results.to_csv(filename, index = True)
        with open(filename, 'wb') as f:
            pickle.dump({'settings': settings, "results": results}, f)
    return results
