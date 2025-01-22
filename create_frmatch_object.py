
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

def create_frmatch_object_adata(adata, cluster_header, nsforest_results): 
    
    # If nsforest_results is str
    if isinstance(nsforest_results, str): 
        nsforest_results = pd.read_csv(nsforest_results)
        
    # Quick NS-Forest results QC
    nsforest_results = nsforest_results.rename(columns = {"NSForest_markers": "markers"})
    if not isinstance(list(nsforest_results["markers"])[0], list): 
        nsforest_results["markers"] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results["markers"]] # Converting str to list format
        if "binary_genes" in list(nsforest_results.columns): 
            nsforest_results["binary_genes"] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results["binary_genes"]] 
    
    # If adata is str
    if isinstance(adata, str): 
        adata = sc.read_h5ad(adata)
    
    # Adding nsforest_results to adata.uns
    dictionary = {cluster_header: {}}
    for index, row in nsforest_results.iterrows(): 
        dictionary[cluster_header][row["clusterName"]] = dict(row[1:])
    adata.uns["nsforest_results"] = dictionary

    return adata

def create_frmatch_object_mtx(cell_by_gene, cluster_header, cluster_labels, nsforest_results, taxonomy = None, save = False): 
    
    adata = ad.AnnData(cell_by_gene)
    adata.obs = cluster_labels.copy()
    adata.obs = adata.obs.astype("category")
    
    adata = create_frmatch_object_adata(adata, cluster_header, nsforest_results)
    
    if taxonomy: 
        dictionary = {"categories_ordered": list(taxonomy)}
        adata.uns["dendrogram_" + cluster_header] = dictionary
        adata.uns["dendrogram_" + cluster_header]["categories_ordered"]

    if save: adata.write_h5ad(save)
    
    return adata

# Converting adata.uns["nsforest_results"]["cluster"] to dataframe
def uns_to_df(dictionary): 
    df = pd.DataFrame(dictionary).transpose().iloc[:, ::-1]
    df["clusterName"] = df.index
    df = df.reset_index(drop = True).iloc[:, ::-1]
    return df
