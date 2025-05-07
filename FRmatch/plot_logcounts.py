
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

def plot_logcounts(adata_query, adata_ref = None): 
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (8, 3))
    a = ax1.hist(adata_query.to_df().unstack(), bins = 20)
    a = ax1.set_title("Histogram of logcounts(adata_query)")
    a = ax1.set_xlabel("logcounts(adata_query)")
    a = ax1.set_ylabel("Frequency")
    a = ax2.hist(adata_ref.to_df().unstack(), bins = 20)
    a = ax2.set_title("Histogram of logcounts(adata_ref)")
    a = ax2.set_xlabel("logcounts(adata_ref)")
    a = ax2.set_ylabel("Frequency")
    return
