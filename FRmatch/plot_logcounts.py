
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib

def plot_logcounts(adata_query, adata_ref = None): 
    """\
    Plots logcounts for input AnnData objects.

    Parameters
    ----------
        adata_query: AnnData
            Annotated data matrix.
        adata_ref: AnnData (default: None)
            Annotated data matrix.
    """
    if adata_ref: 
        fig, axes = plt.subplots(1, 2, figsize = (8, 3))
        p = axes[0].hist(adata_query.to_df().unstack(), bins = 20)
        axes[0].set_title("Histogram of logcounts(adata_query)")
        axes[0].set_ylabel("Frequency")
        axes[0].set_xlabel("logcounts(adata_query)")

        p = axes[1].hist(adata_ref.to_df().unstack(), bins = 20)
        axes[1].set_title("Histogram of logcounts(adata_query)")
        axes[1].set_ylabel("Frequency")
        axes[1].set_xlabel("logcounts(adata_query)")
        
        fig.tight_layout()
        return 
    else: 
        fig, ax = plt.subplots(figsize = (width, height))
        p = ax.hist(adata_query.to_df().unstack(), bins = 20)
        ax.set_title("Histogram of logcounts(adata_query)")
        ax.set_xlabel("logcounts(adata_query)")
        ax.set_ylabel("Frequency")
        return 
