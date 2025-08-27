
# from IPython.core.debugger import set_trace
from .FRtest import FRtest
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib
import random
from random import sample 

def FRtest_subsamp(samp1, samp2, use_cosine = False, subsamp_size = 20, subsamp_iter = 1000, subsamp_seed = False, return_all = False): 
    """\
    Cluster2cluster subsampling samp1 and samp2 to `subsamp_size` `subsamp_iter` times. Returns stat and p_values.

    Parameters
    ----------
        samp1: pd.DataFrame()
            Dataframe representing the query. rows = cells. columns = genes.
        samp2: pd.DataFrame()
            Dataframe representing the reference. rows = cells. columns = genes.
        use_cosine: bool (default: False)
            Whether to use cosine (or Euclidean) distance for tree construction.
        subsamp_size: int (default: 20)
            Number of randomly selected cells per dataset to run tree construction.
        subsamp_iter: int (default: 1000)
            Number of iterations.
        subsamp_seed: bool (default: False)
            Seed for random number generator.
        return_all: bool (default: False)
            Whether to return all FRtest results or just the the median.
    
    Returns
    -------
    out_all_sort: pd.DataFrame()
        Median row after sorted by p-value.
    """
    
    xx = samp1.copy()
    yy = samp2.copy()
    
    m = xx.shape[0]
    n = yy.shape[0]
    
    if subsamp_seed: 
        random.seed(subsamp_seed)

    out_all = pd.DataFrame()
    # for each iteration, run FRtest on a random subset of query and reference
    for b in range(subsamp_iter): 
        mm = sample(range(m), min(subsamp_size, m))
        nn = sample(range(n), min(subsamp_size, n))
        mm = [list(xx.index)[val] for val in mm]
        nn = [list(yy.index)[val] for val in nn]
        xx_B = xx.loc[mm,:]
        yy_B = yy.loc[nn,:]
        out_B = FRtest(xx_B, yy_B, use_cosine = use_cosine)
        out_all = pd.concat([out_all, pd.DataFrame(out_B, index = [0])])
    out_all_sort = out_all.sort_values("p_value", ascending = True)
    
    if return_all: 
        # returning all values
        return out_all_sort 
    else: 
        # returning median row
        return out_all_sort.iloc[int(subsamp_iter/2):int(subsamp_iter/2) + 1,:] 

def FRtest_cell2cluster(samp1, samp2, use_cosine = False, subsamp_size = 20, subsamp_iter = 1000, subsamp_seed = False): 
    """\
    Cell2cluster subsampling samp1 and samp2 to `subsamp_size` `subsamp_iter` times. Keeps the smallest p-value per sample, representing the best FR-match. Returns stat and p_values.

    Parameters
    ----------
        samp1: pd.DataFrame()
            Dataframe representing the query. rows = cells. columns = genes.
        samp2: pd.DataFrame()
            Dataframe representing the reference. rows = cells. columns = genes.
        use_cosine: bool (default: False)
            Whether to use cosine (or Euclidean) distance for tree construction.
        subsamp_size: int (default: 20)
            Number of randomly selected cells per dataset to run tree construction.
        subsamp_iter: int (default: 1000)
            Number of iterations.
        subsamp_seed: bool (default: False)
            Seed for random number generator.
    
    Returns
    -------
    out_max: pd.DataFrame()
        All p-values for each sample and reference cluster pair.
    """

    xx = samp1.copy()
    yy = samp2.copy()
    
    m = xx.shape[0]
    n = yy.shape[0]
    
    if subsamp_seed: 
        random.seed(subsamp_seed)
    
    out_max = {} 
    # for each iteration, run FRtest on a random subset of query and reference
    for b in range(subsamp_iter): 
        mm = sample(range(m), min(subsamp_size, m))
        nn = sample(range(n), min(subsamp_size, n))
        mm = [list(xx.index)[val] for val in mm]
        nn = [list(yy.index)[val] for val in nn]
        xx_B = xx.loc[mm,:]
        yy_B = yy.loc[nn,:]
        out_B = FRtest(xx_B, yy_B, use_cosine = use_cosine)
        # only keeping the smallest p_value for each sample
        for cell in mm: 
            if cell not in out_max: 
                out_max[cell] = out_B
            elif out_max[cell]["p_value"] < out_B["p_value"]: 
                out_max[cell] = out_B
    # reordering to be in the same index (sample) order as samp1
    out_max = pd.DataFrame(out_max).transpose().iloc[:, ::-1].reset_index()
    out_max["index"] = out_max["index"].astype("category")
    out_max["index"] = out_max["index"].cat.set_categories(list(xx.index))
    out_max = out_max.sort_values("index").reset_index(drop = True)
    return out_max
