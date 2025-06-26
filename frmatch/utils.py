
import numpy as np
import pandas as pd
import scanpy as sc

def get_markers(adata, cluster_header, marker_col): 
    """\
    Getting NS-Forest markers in dendrogram order. Removing duplicates. 
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

def compare(df1, df2): 
    df1 = pd.DataFrame(df1)
    df2 = pd.DataFrame(df2)
    df1.index = range(df1.shape[0])
    df1.columns = range(df1.shape[1])
    df2.index = range(df2.shape[0])
    df2.columns = range(df2.shape[1])
    compare = (df1 == df2).copy()
    for index, row in compare.iterrows(): 
        for val in row: 
            if val == False: 
                print("df1", list(df1.iloc[index,:]))
                print("df2", list(df2.iloc[index,:]))
                continue
    print("If nothing else is printed, then dataframes are the same")
