
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

def create_frmatch_object_adata(adata, nsforest_results, marker_col = "binary_genes", additional_markers = [], save = False): 
    # marker_col inputs: NSForest_markers, markers, binary_genes (default: binary_genes)
    
    # Checking adata is ad.AnnData
    if not isinstance(adata, ad.AnnData): 
        print("Error: please pass adata as an anndata object")
        return
    
    # Checking if dendrogram order is provided, if not run scanpy.dendrogram
    cluster_header = list(nsforest_results["cluster_header"])[0]
    if "dendrogram_" + cluster_header not in adata.uns or "categories_ordered" not in adata.uns["dendrogram_" + cluster_header] or len(adata.uns["dendrogram_" + cluster_header]["categories_ordered"]) == 0: 
        sc.tl.dendrogram(adata, cluster_header)
        print(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
    
    # If nsforest_results is str
    if isinstance(nsforest_results, str): 
        if ".csv" in nsforest_results: nsforest_results = pd.read_csv(nsforest_results)
        elif ".pkl" in nsforest_results: nsforest_results = pd.read_pickle(nsforest_results)
    
    # Quick NS-Forest results QC
    if "clusterName" not in list(nsforest_results.columns): 
        print("Error: 'clusterName' not found in nsforest_results")
        return
    elif marker_col not in list(nsforest_results.columns): 
        print("Error: input {marker_col} not found in nsforest_results")
        return
    # converting all values in `marker_col` to lists
    elif not isinstance(list(nsforest_results[marker_col])[0], list): 
        nsforest_results[marker_col] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results[marker_col]] # Converting str to list format
    # one row per clusterName
    if nsforest_results.shape[0] != len(np.unique(nsforest_results["clusterName"])): 
        values = []
        nsforest_new = pd.DataFrame(columns = ["clusterName", marker_col])
        nsforest_new["clusterName"] = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
        for cluster in nsforest_new["clusterName"]: 
            subset = nsforest_results[nsforest_results["clusterName"] == cluster]
            value = []
            for val in subset[marker_col]: 
                value.extend(val)
            values.append(value)
        nsforest_new[marker_col] = values
        nsforest_results = nsforest_new.copy()
        
    # Getting `marker_col` and subsetting adata
    genes = []
    for val in nsforest_results[marker_col]: 
        genes.extend(val)
    
    # Adding additional_markers
    additional_markers = list(set(additional_markers).intersection(set(adata.var_names)))
    genes.extend(additional_markers)
    genes = list(np.unique(genes))
    adata = adata[:,genes]
    
    # If `markers` in nsforest_results, store as ordered list in adata.uns
    col = None
    # from NSForest's nsforesting module
    if "NSForest_markers" in list(nsforest_results.columns): col = "NSForest_markers"
    # from NSForest's evaluating module
    elif "markers" in list(nsforest_results.columns): col = "markers"
    if col: 
        if not isinstance(list(nsforest_results[col])[0], list): 
            nsforest_results[col] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results[col]]
        markers = []
        for val in nsforest_results[col]: 
            val = list(set(val).intersection(set(adata.var.index)))
            val = list(set(val) - set(markers))
            markers.extend(val)
        # saving marker list (ordered according to dendrogram) as adata.uns["markers"]
        adata.uns["markers"] = list(markers)
        adata.uns["markers_per_cluster"] = dict(zip(nsforest_results["clusterName"], nsforest_results[col]))
        
    if save: 
        if isinstance(save, bool): 
            save = f"adata_{marker_col}.h5ad"
        elif ".h5ad" not in save: 
            save = f"{save}.h5ad"
        print(f"Saving anndata object as {save}...")
        adata.write_h5ad(save)

    return adata

def create_frmatch_object_mtx(cell_by_gene, cluster_labels, nsforest_results, marker_col, taxonomy = None, additional_markers = [], save = False): 
    # nsforest_results should have these columns minimum: clusterName, 
    adata = ad.AnnData(cell_by_gene)
    adata.obs = cluster_labels.copy()
    adata.obs = adata.obs.astype("category")
    cluster_header = list(nsforest_results["cluster_header"])[0]
    
    if taxonomy: 
        dictionary = {"categories_ordered": list(taxonomy)}
        adata.uns["dendrogram_" + cluster_header] = dictionary
        adata.uns["dendrogram_" + cluster_header]["categories_ordered"]
    else: 
        adata.uns["dendrogram_" + cluster_header] = {"categories_ordered": []}
    
    adata = create_frmatch_object_adata(adata, nsforest_results, marker_col, additional_markers = additional_markers, save = save)
    return adata

# Converting adata.uns["nsforest_results"]["cluster"] to dataframe
def uns_to_df(dictionary): 
    df = pd.DataFrame(dictionary).transpose().iloc[:, ::-1]
    df["clusterName"] = df.index
    df = df.reset_index(drop = True).iloc[:, ::-1]
    return df
