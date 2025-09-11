
'''
srun --partition=ind-shared --pty --account=jcl125 --nodes=1 --ntasks-per-node=8 --mem=243G -t 8:00:00 --wait=0 --export=ALL /bin/bash

conda activate environment
cd frmatch
python3
'''

import os
import sys
sys.path.insert(0, os.path.abspath("./FRmatch"))
import FRmatch
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import time

# This is the Python code for running the FRmatch pipeline cell2cluster matching CellRef clusters to HLCA clusters
# Files required to run this script: HLCA NSForest results, CellRef NSForest results

# Defining additional markers: HLCA and CellRef binary_genes union
additional_markers = []
nsforest_results = pd.read_csv("../nsforest/NSForest/outputs_hlca_core_nsforest_markers/hlcacore_markers_ann_fine_level.csv")
nsforest_results["binary_genes"] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results["binary_genes"]]
for val in nsforest_results["binary_genes"]: 
    additional_markers.extend(val)

nsforest_results = pd.read_csv("../nsforest/NSForest/outputs_cellref_nsforest_markers/celltype_level3_results.csv")
nsforest_results["binary_genes"] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results["binary_genes"]]
for val in nsforest_results["binary_genes"]: 
    additional_markers.extend(val)

additional_markers = list(np.unique(additional_markers))

type_ = "cell2cluster"
# bin_ = 0 # bin number 0-5
bin_ = False
fold = 10 
output_folder = "tutorial_hlca_cellref"
output_filename = f"results_CellReftoHLCA_{type_}_fold{fold}"
# output_filename = f"results_HLCAtoCellRef_{type_}_bin{bin_}_fold{fold}"

# Preparing HLCA anndata object
adata_hlca = sc.read_h5ad("../nsforest/NSForest/demo_data/hlca_core.h5ad")
adata_hlca.var = adata_hlca.var.reset_index()
adata_hlca.var.index = adata_hlca.var["feature_name"]
cluster_header_hlca = "ann_finest_level"
nsforest_results = pd.read_csv("../nsforest/NSForest/outputs_hlca_core_nsforest_markers/hlcacore_markers_ann_fine_level.csv")
nsforest_results["binary_genes"] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results["binary_genes"]]
nsforest_results["cluster_header"] = cluster_header_hlca
save = "data/adata_hlca_frmatch.h5ad"
adata_hlca = FRmatch.create_frmatch_object_adata(adata_hlca, nsforest_results, marker_col = "binary_genes", additional_markers = additional_markers, save = False)

# Preparing CellRef anndata object
adata_cellref = sc.read_h5ad("../nsforest/NSForest/demo_data/LungMAP_HumanLung_CellRef.v1.1.h5ad") 
adata_cellref.var.index = adata_cellref.var["features"]
cluster_header_cellref = "celltype_level3_fullname"
adata_cellref.obs[cluster_header_cellref] = adata_cellref.obs[cluster_header_cellref].astype('category')
dictionary = dict(zip(adata_cellref.obs["celltype_level3"], adata_cellref.obs[cluster_header_cellref]))
nsforest_results = pd.read_csv("../nsforest/NSForest/outputs_cellref_nsforest_markers/celltype_level3_results.csv")
nsforest_results["binary_genes"] = [val.replace("'", "").replace("[", "").replace("]", "").split(", ") for val in nsforest_results["binary_genes"]]
nsforest_results["clusterName"] = [dictionary[val] for val in nsforest_results["clusterName"]]
nsforest_results["cluster_header"] = cluster_header_cellref
save = "data/adata_cellref_frmatch.h5ad"
adata_cellref = FRmatch.create_frmatch_object_adata(adata_cellref, nsforest_results, marker_col = "binary_genes", additional_markers = additional_markers, save = False)

for col in adata_hlca.obs.columns: 
    if "ann_finest_level" not in col: 
        del adata_hlca.obs[col]

for col in adata_cellref.obs.columns: 
    if "celltype_level3_fullname" not in col: 
        del adata_cellref.obs[col]

print(adata_hlca)
print(dict(adata_hlca.obs["ann_finest_level"].value_counts()))
print(adata_cellref)
print(dict(adata_cellref.obs["celltype_level3_fullname"].value_counts()))

# Binning clusters by size
# Returns a list of length {partitions} with the clusters divided
def binning(adata, cluster_header, partitions = 6): 
    # Getting value counts sorted by size
    dictionary = dict(adata.obs[cluster_header].value_counts())
    # Getting minimum number of cells for bin to close
    minimum = int(sum(dictionary.values())/partitions)
    bins = []
    values, count = [], 0
    for key in dictionary.keys(): 
        # print(key, dictionary[key])
        # If the current count > minimum, add values to bin and reset
        if count > minimum: 
            bins.append(values)
            values, count = [], 0
        values.append(key)
        count += dictionary[key]
    bins.append(values)
    return bins

# Note: binning is not required but recommended because the HLCA and CellRef AnnData are large files

# Separating cellref into folds
folds = pd.read_csv("cellref_tenfold_idx 2.csv") # cols: fold_{1-10} 0-indexing
## take care of small clusters, except Chondrocyte
ind_small, index = [], 0
for val in adata_cellref.obs[cluster_header_cellref]: 
    if val in ["Megakaryocyte/Platelet", "Innate lymphoid cell"]: 
    # if val in ["Megakaryocyte/Platelet", "Innate lymphoid cell", "Chondrocyte"]: 
        ind_small.append(index)
    index += 1
folds = list(folds[f"fold_{fold}"])
# Adding ind_small to fold 1, removing from all other folds
if fold == 1: 
    folds.extend(ind_small)
    folds = list(set(folds))
else: 
    folds = list(set(set(folds) - set(ind_small)))
adata_cellref[folds,:].obs[cluster_header_cellref].value_counts()
adata_cellref = adata_cellref[folds,:]

print(adata_hlca)
print(adata_cellref)

# Normalizing
adata_hlca = FRmatch.normalization(adata_hlca, "ann_finest_level", norm_by = "mean")
adata_cellref = FRmatch.normalization(adata_cellref, "celltype_level3_fullname", norm_by = "mean")

print(adata_hlca)
print(adata_cellref)

results_CellReftoHLCA = FRmatch.FRmatch_cell2cluster(adata_cellref, adata_hlca, cluster_header_query = "celltype_level3_fullname", cluster_header_ref = "ann_finest_level", use_cosine = True, filter_size = 10, subsamp_size = 10, subsamp_iter = 1000, verbose = 1, save = f"{output_folder}/{output_filename}")
