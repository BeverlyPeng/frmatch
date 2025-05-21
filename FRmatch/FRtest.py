
from IPython.core.debugger import set_trace
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def FRtest(samp1, samp2, use_cosine = False, 
           plot_mst = False, colors = ["#F0E442", "#56B4E9"], label_names = ["Sample 1", "Sample 2"], 
           vertex_size = 50, edge_width = 0.5, title = "Minimum spanning tree"): 
    """\
    Generating minimum spanning tree (MST) of 2 samples. 
    """
    xx = samp1
    yy = samp2
    
    ## get sample sizes
    m = xx.shape[0]
    n = yy.shape[0]
    N = m + n
    
    ##--- OPTION 2: BINARY VALUES ---##
    
    ## column-wise combined data matrix
    vectors = np.concatenate((xx, yy), axis = 0)
    
    ## euclidean distance
    G = nx.Graph()
    for i in range(len(vectors)):
        if sum(vectors[i]) == 0: 
            # Add warning
            vectors[i] = np.array([1e-10]*len(vectors[i]))
        for j in range(i + 1, len(vectors)):
            if sum(vectors[j]) == 0: 
                vectors[j] = np.array([1e-10]*len(vectors[j]))
            # Calculate the Euclidean distance between vectors[i] and vectors[j]
            if not use_cosine: 
                distance = np.linalg.norm(np.array(vectors[i]) - np.array(vectors[j]))
            ##--- OPTION 1: CONSINE DISTANCE ---##
            else: 
                distance = 1 - np.dot(np.array(vectors[i]), np.array(vectors[j]))/(np.linalg.norm(np.array(vectors[i]))*np.linalg.norm(np.array(vectors[j])))
            G.add_edge(i, j, weight=distance)
    
    # MST
    mst = nx.minimum_spanning_tree(G, algorithm='kruskal')
    # Convert MST to binary values in numpy
    myMST = np.ceil(nx.adjacency_matrix(mst).toarray())
    myMST[myMST != 0] = 1
    myMST = myMST.astype(int)
    
    ## count total runs (i.e. total subgraphs)
    bottomleft = myMST[m:,:m] # myMST[(m+1):N,1:m]
    runs = sum(sum(bottomleft)) + 1
    
    ## count subgraphs for each sample
    bottomright = myMST[m:,m:] # myMST[(m+1):N,(m+1):N]
    G = nx.Graph(bottomright)
    runs_samp2 = nx.number_connected_components(G)
    topleft = myMST[:m,:m] # myMST[1:m,1:m]
    G = nx.Graph(topleft)
    runs_samp1 = nx.number_connected_components(G)
    
#     if runs != runs_samp1 + runs_samp2: 
        

    ## calculate common nodes
    xsum = myMST.sum(axis = 0)
    C = sum(xsum*(xsum-1))/2
    ## calculate mean and variance
    mu = 2*m*n/N + 1
    sigma_sq = (2*m*n/(N*(N-1)))*((2*m*n-N)/N+(C-N+2)*(N*(N-1)-4*m*n+2)/((N-2)*(N-3)))
    ## the standardized FR-stat and p-value
    stat = (runs-mu)/np.sqrt(sigma_sq)
    p_value = norm.cdf(stat)
    
#     if runs in [16, 17]: 
#         print("C", C)
#         print("runs", runs, runs_samp1, runs_samp2)
#         print(myMST.shape, topleft.shape, bottomright.shape)
#         print("m, n, N", m, n, N)
#         print("mu", mu)
#         print("sigma_sq", sigma_sq)
#         print("stat", stat)
#         print("p_value", p_value)
#         print()

    ## plot
    if plot_mst: 
        plt.figure()
        # https://networkx.org/documentation/stable/reference/generated/networkx.drawing.nx_pylab.draw_networkx.html
        labels = {}
        for i in range(m): labels[i] = colors[0]
        for i in range(m, N): labels[i] = colors[1]
    #     pos = nx.bfs_layout(mst, 20)
    #     nx.draw_networkx(mst, pos = pos, node_color = labels.values(), with_labels = True, node_size = 50, font_size = 10, label = label_names[0])
        nx.draw_networkx(mst, with_labels = False, node_size = vertex_size, node_color = labels.values(), 
                         width = edge_width, edge_color = "gray", font_size = 10, #label = "test", #label_names[0] 
                         edgecolors = "black", linewidths = 0.5, 
                        )
        plt.title(title)
        dict_elements = {"marker": 'o', "linewidth": 0, "markeredgewidth": 0.5, "color": 'black', "markersize": 7}
        legend_elements = [Line2D([0], [0], label = f"{label_names[0]} ({m})", markerfacecolor=colors[0], **dict_elements), 
                           Line2D([0], [0], label = f"{label_names[1]} ({n})", markerfacecolor=colors[1], **dict_elements)]
        plt.legend(handles = legend_elements)
        plt.show()
    
    # output
    return {"runs": runs, "runs_samp1": runs_samp1, "runs_samp2": runs_samp2, "stat": stat, "p_value": p_value}
