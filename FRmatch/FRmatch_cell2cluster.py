
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

def FRmatch_cell2cluster(query, ref, cluster_header_query, cluster_header_ref, feature_selection = "reference_markers", 
                         use_cosine = True, filter_size = 5, subsamp_size = 10, subsamp_iter = 2000, subsamp_seed = 916, 
                         subsamp_iter_custom = False, subsamp_iter_custom_k = 5, 
                         prefix = ["query", "ref"], verbose = 0, save = False): 
    """\
    FRmatch cell to cluster matching.

    Parameters
    ----------
        query: AnnData
            Query object.
        ref: AnnData
            Reference object.
        cluster_header_query: str
            Column in `query.obs` storing cell annotation.
        cluster_header_ref: str
            Column in `ref.obs` storing cell annotation.
        feature_selection: str (default: "reference_markers")
            Method to select genes
        use_cosine: bool (default: True)
            Whether to use cosine distance vs Euclidean distance for tree construction.
        filter_size: int (default: 5)
            Minimum cluster size.
        subsamp_size: int (default: 10)
            Number of cells per dataset to run tree construction.
        subsamp_iter: int (default: 1000)
            Number of iterations
        subsamp_seed: bool (default: False)
            Seed for random number generator.
        subsamp_iter_custom: bool (default: False)
            Whether to use custom subsamp_iter to ensure all cells are adequately sampled.
        subsamp_iter_custom_k: int (default: 5)
            Number of times * subsamp_iter.
        prefix: list (default: ["query", "ref"])
            Labels for query and reference data.
        verbose: int (default: 0)
            Controls print statements.
        return_all: bool (default: False)
            Whether to return all FRtest results or the median.
        save: bool | str (default: False)
            Whether to save results as pkl. If string, save as string.
    
    Returns
    -------
    dictionary
        FRmatch results with keys ["settings", "p_values", "stat"]
    """
    # Saving settings as dictionary
    settings = {"query": prefix[0], "ref": prefix[1], "cluster_header_query": cluster_header_query, "cluster_header_ref": cluster_header_ref, "feature_selection": feature_selection, "use_cosine": use_cosine, "filter_size": filter_size, "subsamp_size": subsamp_size, "subsamp_iter": subsamp_iter, "subsamp_seed": subsamp_seed, "subsamp_iter_custom": subsamp_iter_custom, "subsamp_iter_custom_k": subsamp_iter_custom_k, "save": save}
    
    if save: 
        if not isinstance(save, str): 
            filename = f"frmatch_results_cell2cluster_{prefix[0]}to{prefix[1]}.pkl"
        elif ".pkl" not in str(save): 
            filename = save + ".pkl"
        else: 
            filename = save
    
    ## feature selection
    ## use query gene space
    if feature_selection == "query_genes": 
        marker_genes = list(query.var_names)
        marker_genes_common = list(set(marker_genes) & set(ref.var_names))
        if verbose > 0: 
            print(f"Feature selection: marker genes of the reference experiment.")
            print(f"\t{len(marker_genes_common)} out of {len(marker_genes)} query genes are common in the reference experiment.")
        if len(marker_genes) != len(marker_genes_common) and verbose > 0: 
            print(f"{sorted(list(set(marker_genes) - set(marker_genes_common)))} not found.")
    ## use reference markers
    elif feature_selection == "reference_markers": 
        marker_genes = ref.uns["markers"]
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
    
    ## filtering small clusters
    if verbose > 0: 
        print(f"Filtering small clusters: query and reference clusters with less than {filter_size} cells are not considered.")
        
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
    
    ## extract cluster info from adata.obs
    membership_query = query_reduced.obs[cluster_header_query]
    membership_ref = ref_reduced.obs[cluster_header_ref]
    
    ncluster_query = len(np.unique(membership_query))
    ncluster_ref = len(np.unique(membership_ref))
    if verbose > 0: print(f"Comparing {ncluster_query} query clusters with {ncluster_ref} reference clusters.")
    
    ## FR comparison between two experiments
    results = pd.DataFrame()
    for query_cluster in tqdm(np.unique(membership_query), desc = "FRmatch"): 
        warning_message = True

        for ref_cluster in np.unique(membership_ref): 
            if verbose > 0: print("QUERY: ", query_cluster, "ref: ", ref_cluster)
            
            # subsetting query and ref to query_cluster and ref_cluster
            query_df = query_X_filt[query_X_filt["cluster"] == query_cluster]
            del query_df["cluster"]
            ref_df = ref_X_filt[ref_X_filt["cluster"] == ref_cluster]
            del ref_df["cluster"]
            if verbose > 0: print(query_df.shape, ref_df.shape)
            
            if subsamp_iter_custom: 
                subsamp_iter = max(subsamp_iter, subsamp_iter_custom_k * query_df.shape[0])
                
            df = FRmatch.FRtest_cell2cluster(query_df, ref_df, subsamp_size = subsamp_size, subsamp_iter = subsamp_iter, subsamp_seed = subsamp_seed)
            
            if df.shape[0] != query_df.shape[0] and warning_message: 
                print(f"warning: Not all cells are randomly sampled. Consider increasing the number of iterations specified in the 'subsamp_iter' argument. query_cluster: {query_cluster} ref_cluster: {ref_cluster}")
                warning_message = False
            
            df["query_cluster"] = query_cluster
            df["ref_cluster"] = ref_cluster
            results = pd.concat([results, df])
            if save: 
                with open(filename, 'wb') as f:
                    pickle.dump({'settings': settings, "results": results}, f)
            
    if save: 
        with open(filename, 'wb') as f:
            pickle.dump({'settings': settings, "results": results}, f)
    return results
