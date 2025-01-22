
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

# the next step is to set up the iterative subsampling based on FRtest
def FRtest_subsamp(samp1, samp2, subsamp_size = 20, subsamp_iter = 1000, subsamp_seed = 916, return_all = False): 
    """\
    Subsampling samples to `subsamp_size` `subsamp_iter` times. Returns dataframe of stat and p_values.  
    """
    ## data input matrices: rows = cells, columns = genes
    xx = samp1.copy()
    yy = samp2.copy()
    
    ## get original sample sizes
    m = xx.shape[0]
    n = yy.shape[0]
    
    # setting seed 
#     print("before", random.random())
    random.seed(subsamp_seed)
#     print("after", random.random())
#     random.seed(subsamp_seed)
#     print("after", random.random())
    out_all = pd.DataFrame()
    for b in range(subsamp_iter): 
#         print("during", random.random())
        mm = sample(range(m), min(subsamp_size, m))
        nn = sample(range(n), min(subsamp_size, n))
        mm = [list(xx.index)[val] for val in mm]
        nn = [list(yy.index)[val] for val in nn]
        xx_B = xx.loc[mm,:]
        yy_B = yy.loc[nn,:]
        out_B = FRtest(xx_B, yy_B)
        out_all = pd.concat([out_all, pd.DataFrame(out_B, index = [0])])
    out_all_sort = out_all.sort_values("p_value", ascending = True)
    
    if return_all: 
        return out_all_sort
    else: 
        return out_all_sort.iloc[int(subsamp_iter/2):int(subsamp_iter/2) + 1,:] # returning median row

def FRtest_cell2cluster(samp1, samp2, subsamp_size = 20, subsamp_iter = 1000, subsamp_seed = 916): 
    ## data input matrices: rows = cells, columns = genes
    xx = samp1.copy()
    yy = samp2.copy()
    
    ## get original sample sizes
    m = xx.shape[0]
    n = yy.shape[0]
    
    # setting seed 
#     print("before", random.random())
    random.seed(subsamp_seed)
    
    out_max = {} # E18_1_Nuclei_NeuNP_H200_1030_MTG_layer1_BCH7: {'runs': 2, 'runs_samp1': 1, 'runs_samp2': 1, 'stat': -4.307588553812465, 'p_value': 8.25220137304871e-06}
    for b in range(subsamp_iter): 
#         print("during", random.random())
        mm = sample(range(m), min(subsamp_size, m))
        nn = sample(range(n), min(subsamp_size, n))
        mm = [list(xx.index)[val] for val in mm]
        nn = [list(yy.index)[val] for val in nn]
        xx_B = xx.loc[mm,:]
        yy_B = yy.loc[nn,:]
        out_B = FRtest(xx_B, yy_B)
        for cell in mm: 
            if cell not in out_max: 
                out_max[cell] = out_B
            elif out_max[cell]["p_value"] < out_B["p_value"]: 
                out_max[cell] = out_B
        # https://www.datacamp.com/tutorial/pipe-r-tutorial
        # out.cell2cluster[mm] %<>% pmax(out.B["p.value"], na.rm=TRUE) 
        # pmax: parallel maxima, meaning replace the original value in out.cell2cluster if out.B["p.value"] is larger
#         out_all = pd.concat([out_all, pd.DataFrame(out_B, index = [0])])
#     out_all_sort = out_all.sort_values("p_value", ascending = True)
    out_max = pd.DataFrame(out_max).transpose().iloc[:, ::-1].reset_index()
    out_max["index"] = out_max["index"].astype("category")
    out_max["index"] = out_max["index"].cat.set_categories(list(xx.index))
    out_max = out_max.sort_values("index").reset_index(drop = True)
    return out_max