
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

def FRmatch_cluster2cluster(query, ref, cluster_header_query, cluster_header_ref, marker_genes = [], 
                            use_cosine = True, filter_size = 5, subsamp_size = 20, subsamp_iter = 1000, subsamp_seed = False, 
                            prefix = ["query", "ref"], verbose = 0, return_all = False, save = False): 
    """\
    Creating an AnnData object to run FRmatch on.

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
        marker_genes: list (default: [])
            Genes to subset instead of ref.uns["markers"].
        use_cosine: bool (default: True)
            Whether to use cosine distance vs Euclidean distance for tree construction.
        filter_size: int (default: 5)
            Minimum cluster size.
        subsamp_size: int (default: 20)
            Number of cells per dataset to run tree construction.
        subsamp_iter: int (default: 1000)
            Number of iterations
        subsamp_seed: bool (default: False)
            Seed for random number generator.
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
    settings = {"query": prefix[0], "ref": prefix[1], "cluster_header_query": cluster_header_query, "cluster_header_ref": cluster_header_ref, "marker_genes": marker_genes, "use_cosine": use_cosine, "filter_size": filter_size, "subsamp_size": subsamp_size, "subsamp_iter": subsamp_iter, "subsamp_seed": subsamp_seed, "return_all": return_all, "save": save}
    
    if save: 
        if not isinstance(save, str): 
            filename = f"frmatch_results_cluster2cluster_{prefix[0]}to{prefix[1]}.pkl"
        elif ".pkl" not in str(save): 
            filename = save + ".pkl"
        else: 
            filename = save

    ## find common marker genes from ref
    if len(marker_genes) == 0: 
        marker_genes = ref.uns["markers"]
    
    marker_genes_common = []
    for marker in marker_genes: 
        if marker in list(query.var_names): 
            marker_genes_common.append(marker)
    if verbose > 0: 
        print(f"Feature selection: {len(marker_genes_common)} out of {len(marker_genes)} reference marker genes are presented in the query experiment.")
    if len(marker_genes) != len(marker_genes_common) and verbose > 0: 
        print(f"{sorted(list(set(marker_genes) - set(marker_genes_common)))} not found.")
    
    ## filtering small or low fscore clusters
    if verbose > 0: 
        print(f"Filtering small clusters: query and reference clusters with less than {filter_size} cells are not considered.")
    query_shape = query.shape
    ref_shape = ref.shape
    query = FRmatch.filter_cluster(query, cluster_header_query, filter_size) #only filter on size, not fscore
    ref = FRmatch.filter_cluster(ref, cluster_header_ref, filter_size)
    if verbose > 0: 
        print(f"filtered query from {query_shape} -> {query.shape}")
        print(f"filtered ref from {ref_shape} -> {ref.shape}")
    
    ## get expr data and dimension reduction by selecting common marker genes
    query_X_filt = query[:, marker_genes_common].to_df()
    query_X_filt["cluster"] = query.obs[cluster_header_query]
    ref_X_filt = ref[:, marker_genes_common].to_df()
    ref_X_filt["cluster"] = ref.obs[cluster_header_ref]
    
    ## extract cluster info from sce.objects
    query_clusters = list(query.obs[cluster_header_query])
    ref_clusters = list(ref.obs[cluster_header_ref])

    ## FR comparison between two experiments
    ## prepare data for each combination of query cluster and ref cluster pair
    results = pd.DataFrame()
    for query_cluster in tqdm(np.unique(query_clusters), desc = "FRmatch"): 

        for ref_cluster in np.unique(ref_clusters): 
            if verbose > 0: print("QUERY: ", query_cluster, "ref: ", ref_cluster)

            # subsetting query and ref to query_cluster and ref_cluster
            query_df = query_X_filt[query_X_filt["cluster"] == query_cluster]
            del query_df["cluster"]
            ref_df = ref_X_filt[ref_X_filt["cluster"] == ref_cluster]
            del ref_df["cluster"]
            if verbose > 0: print(query_df.shape, ref_df.shape)

            df = FRmatch.FRtest_subsamp(query_df, ref_df, use_cosine = use_cosine, subsamp_size = subsamp_size, subsamp_iter = subsamp_iter, subsamp_seed = subsamp_seed, return_all = return_all)
            df["query_cluster"] = query_cluster
            df["ref_cluster"] = ref_cluster
            results = pd.concat([results, df])

        if save: 
            with open(filename, 'wb') as f:
                pickle.dump({'settings': settings, "p_values": results.pivot(index = "ref_cluster", columns = "query_cluster", values = "p_value"), "stat": results.pivot(index = "ref_cluster", columns = "query_cluster", values = "stat")}, f)
        
    results_pval = results.pivot(index = "ref_cluster", columns = "query_cluster", values = "p_value")
    if save: 
        print(f"Saving frmatch results as... {filename}")
        with open(filename, 'wb') as f:
            pickle.dump({'settings': settings, "p_values": results.pivot(index = "ref_cluster", columns = "query_cluster", values = "p_value"), "stat": results.pivot(index = "ref_cluster", columns = "query_cluster", values = "stat")}, f)

    return {'settings': settings, "p_values": results.pivot(index = "ref_cluster", columns = "query_cluster", values = "p_value"), "stat": results.pivot(index = "ref_cluster", columns = "query_cluster", values = "stat")}
