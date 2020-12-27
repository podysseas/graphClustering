# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 17:10:24 2020

@author: ody
"""

"""

"""
from random import randint

import sys
import numpy as np
import os
import pandas as pd
import editdistance
import time
import igraph as ig
import seaborn as sns
import statistics
from sklearn.datasets import make_circles
from sklearn.neighbors import kneighbors_graph
from sklearn.cluster import SpectralClustering

def add_weights_as_list(g,distance):
    
    '''
        Input : graph 
                distance/similarity matrix
        -----------------------------
        Output: graph with attribute weights completed
    
    '''
    g = graph_threshold
    

    weights = []

    edges = g.get_edgelist()
    
    for i in range(0,len(edges[:])):
        # find weight from distance matrix
        i_index = edges[i][0]
        j_index = edges[i][1]
        weights.append(distance[i_index,j_index])
        
    g.es['weight'] = weights 
    return g

def all_same(items):
    return all(x == items[0] for x in items)

def Plot_Graph(graph,label):
    #graph.vs["label"] = label
    color_dict = []
    values, counts = np.unique(label, return_counts=True)

    n = len(values)+1 #clusters start from 0 or 1!
    for i in range(n):
        color_dict.append('#%06X' % randint(0, 0xFFFFFF))
    
    np.unique(label)
    graph.vs["color"] = [color_dict[l] for l in label]
    return graph

def Plot_Distribution(distance):
    '''
        Input : similarity matrix
        -------------------------
        Output : figure object
    '''
    A = distance.flatten()
    
    t = sns.distplot(A, hist=True, kde=False,
                     bins=int(300 / 5), color="darkblue",
                     hist_kws={'edgecolor': 'black'},
                     kde_kws={'linewidth': 4})
    figure = t.get_figure()
    return figure

def calculate_distance(data):
    
    '''
        Input  : a matrix with string sequences 
        --------------------------------------
        Output : a similarity matrix normalized to [0,1]
    '''
    limit = len(data)
    distance = np.zeros((len(data), len(data)))

    for i in range(0, len(distance)):
        limit = len(filtered_data)

        for j in range(limit):
            string1 = data.iloc[i, 1]  # take the second column
            string2 = data.iloc[j, 1]
            temp_distance = 1 - editdistance.eval(string1, string2) / len(string1)
            # fast implementation of Levenshtein
            # https://github.com/roy-ht/editdistance
            temp_distance = format(temp_distance, '.3f')
            distance[i, j] = temp_distance
            distance[j, i] = temp_distance

        limit = limit - 1
    return distance


def Choose_Threshold_Graph(distance, threshold):
    
    ''' 
        Input  :    similarity matrix
                    threshold
        -------------------------
        Output :    graph_threshold
    '''
    size = distance.shape[0]
    adj_matrix = np.ones((size, size), dtype=int)

    np.fill_diagonal(adj_matrix, 0)  # reject self - loops

    graph_threshold = ig.Graph.Adjacency((distance > threshold).tolist())  # create a graph from adjacency matrix

    return graph_threshold


def normalized_weighted_graph(distance):
    
    ''' 
        Input : similarity matrix
        -------------------------
        Output: normalized_vertex_graph
    '''
    size = distance.shape[0]
    global adj_matrix_boolean
    adj_matrix_boolean = np.ones((size, size), dtype=int)
    vertex_threshold = np.ones(size)
    for i in range(0, size):
        temp_threshold = statistics.mean(distance[i, :])
        vertex_threshold[i] = format(temp_threshold, '.4f')

        adj_matrix_boolean[i, :] = distance[i, :] > vertex_threshold[i]

    graph_normalized = ig.Graph.Adjacency(adj_matrix_boolean.tolist())  # create a graph from adjacency matrix

    return graph_normalized


def shared_nearest_neighbors(adj_matrix_boolean,threshold_neighbors):
    '''
        Input:  adjacency matrix of the graph 
                threshold for the neighbors
        -------------------------------------
        Output: graph_shared
                label
                
    '''
    
    size = adj_matrix_boolean.shape[0]
    adj_matrix_shared = np.zeros((size, size))
    max_neighbors     = np.zeros(size)    
    label = np.zeros(size,dtype=int)
    count_label = 1
    
    for i in range(0, size):

        if i == 0:
            label[i] = 1
            count_label = count_label + 1  # new label found

        for j in range(0, size):
            # neighbors = adj_matrix_boolean[i,:] == adj_matrix_boolean[j,:]
            temp_i = adj_matrix_boolean[i, :]
            temp_j = adj_matrix_boolean[j, :]
            neighbors = np.array(np.where(temp_i==temp_j ))

            if neighbors.size > threshold_neighbors:

                if label[i] == 0:  # label is empty for i

                    if label[j] != 0:  # label j is not empty

                        # use label j 
                        label[i] = label[j]

                    if label[j] == 0:
                        # not label assigned
                        label[j] = count_label
                        label[i] = count_label
                        count_label = count_label + 1 # new label
                        if neighbors.size > max_neighbors[i]:
                            max_neighbors[i] = neighbors.size
                
                if label[i] !=0:
                    # if i has already a label
                    
                    if  neighbors.size > max_neighbors[i]:
                        # new label 
                        max_neighbors[i] = neighbors.size
                        label[i] = label[j]
                        
                adj_matrix_shared[i, j] = 1
                adj_matrix_shared[j, i] = 1

        
    graph_shared = ig.Graph.Adjacency(adj_matrix_shared.tolist())  # create a graph from adjacency matrix
    
    
    return graph_shared,label;


if __name__ == "__main__":
    
    FILE = '2_IMGT-gapped-nt-sequences.txt'
    COLUMN = 'V.D.J.REGION'
    SEQUENCE_NUMBER = 'Sequence.number'
    
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    data = pd.read_csv(os.path.join(__location__, FILE), sep='\t')
    filtered_data = data[['Sequence.number', COLUMN]].dropna()  # drop nan
    
    ''' this block calculates distance '''
    start = time.time()
    distance = calculate_distance(filtered_data)
    end = time.time()
    elapsed_time = format(end - start, '.3f')
    print('The execution time is {0} seconds\n'.format(elapsed_time))
    
    ''' ***************************    '''
    layout              ="kk"

    #--------------------------------------------------------
    figure = Plot_Distribution(distance)
    figure.savefig('distribution_of_similarities.png', dpi=400)
    #--------------------------------------------------------
    graph_threshold = Choose_Threshold_Graph(distance, threshold=0.86)
    #graph_threshold = Plot_Graph(graph_threshold)
    ig.plot(graph_threshold, "graph_threshold.png", layout=layout, bbox=(500, 500), margin=50, inline='None')

    #--------------------------------------------------------
    graph_normalized = normalized_weighted_graph(distance)
    #graph_normalized = Plot_Graph(graph_normalized)
    ig.plot(graph_normalized, "graph_normalized.png", layout=layout, bbox=(500, 500), margin=50, inline='None')

    #--------------------------------------------------------
    temp_graph          = graph_normalized
    adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
    threshold_neighbors = 300
    graph_shared,label  = shared_nearest_neighbors(adj_matrix_boolean,threshold_neighbors)
    graph_shared        = Plot_Graph(graph_shared,label)
    layout              ="kk"
    ig.plot(graph_shared, "graph_shared.png", layout=layout, bbox=(500, 500), margin=50, inline='None')
    #--------------------------------------------------------
    
    
    '''
        TO DO : 
            - check the algorithms
            
            Walktrap  --- checked
            
            Fastgreedy -- checked
    
    ''' 

    g = graph_threshold
    g = add_weights_as_list(g,distance)

    wtrap = g.community_walktrap(weights = g.es['weight'],steps = 30)
    lab =np.array(wtrap.as_clustering()).tolist()

    #g.to_undirected()
    #fastgreedy = g.community_fastgreedy(weights = g.es['weight'])    
    #lab =np.array(fastgreedy.as_clustering()).tolist()

    
    labels = np.zeros(len(distance),dtype="int")
    n_clusters = len(lab[:])
    for k in range(0,n_clusters):
        
        for l in range(0,len(lab[k][:])) :
            
            # k is the cluster
            # l is the vertex node 
            i_index = lab[k][l]
            labels[i_index] = k
    
    g        = Plot_Graph(g,labels)

    layout              ="kk"

    ig.plot(g, "graph_plot_3.png", layout=layout, bbox=(500, 500), margin=50, inline='None')
    
    
    
    #matrix_to_spectral = graph_normalized.get_adjacency()[:, :]  # int 0 1
    # matrix_to_spectral = graph_threshold.get_adjacency()[:, :]  # int 0 1
    
    # matrix_to_spectral = np.array(matrix_to_spectral.data)  # int
    # ''' calculate the degree of the matrix '''
    # d = np.diag(matrix_to_spectral.sum(axis=1))
    
    # ''' Laplacian Matrix '''
    # l = np.array(d - matrix_to_spectral)
    # eig_L, eig_vec_L = np.linalg.eig(l)
    
    # i = np.where(eig_L < 0.1)
    
    #sc               = SpectralClustering(5, affinity='precomputed', n_init=100,assign_labels='discretize')   
    
    #label           = sc.fit_predict(adjacency_matrix)
    

    #ig.plot(temp_graph, "graph.png", layout=layout, bbox=(500, 500), margin=50, inline='None')
                


