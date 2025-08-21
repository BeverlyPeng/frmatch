
import numpy as np
import pandas as pd
import scanpy as sc

def get_markers(adata, cluster_header, marker_col): 
    """\
    Getting NS-Forest markers in dendrogram order. Removing duplicates. 

    Parameters
    ----------
        adata: AnnData
            Annotated data matrix.
        cluster_header: str
            Column in adata.obs storing cell annotation.
        marker_col: str
            Column in nsforest_results to subset AnnData. Other options include ['NSForest_markers', 'markers', 'binary_genes'].
    """
    markers = []
    cluster_order = []
    # Checking if custom cluster order is available
    if "cluster_order" in list(adata.uns.keys()): 
        cluster_order = list(adata.uns["cluster_order"])
    # Checking if dendrogram available
    elif "dendrogram_" + cluster_header in list(adata.uns.keys()): 
        cluster_order = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
    else: 
        cluster_order = list(adata.uns["nsforest_results"][cluster_header])
    # Only including specified clusters
    for cluster in cluster_order: 
        if cluster not in list(adata.uns["nsforest_results"][cluster_header]): continue
        values = list(adata.uns["nsforest_results"][cluster_header][cluster][marker_col])
        for marker in markers: 
            if marker in values: 
                values.remove(marker)
        markers.extend(values)
    return markers
