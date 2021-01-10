"""
Created on Sat Sep 19 17:10:24 2020

@author: ody
"""

"""

"""
from random import randint
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import os
import pandas as pd
import editdistance
import time
import igraph as ig
import seaborn as sns
import statistics
import datetime


# from sklearn.datasets import make_circles
# from sklearn.neighbors import kneighbors_graph
# from sklearn.cluster import SpectralClustering

    
def Calculate_distance_metric_inside_cluster(distance,clusters):
    
    clusters = np.array(clusters).tolist()
    
    average_cluster = np.zeros(shape =len(clusters))
    variance_cluster = np.zeros(shape =len(clusters))

    for clust in range(len(clusters)):
        
        nodes = clusters[clust]
       
        if len(nodes) > 1:
            
            sum = 0;
            sum_var = 0
            count = 0
            for i in range(len(nodes)):
                
                k = i + 1
                for j in range(k,len(nodes)):
            
                    if nodes[i] !=nodes[j]:
                        
                        temp = distance[nodes[i],nodes[j]]
                        sum  = sum + temp 
                        
                        count = count + 1 
                        #print(nodes[i],nodes[j])
                    
            average_cluster[clust]  = sum/count
            #-------------------- VARIANCE --------
            sum_var = 0
            count = 0
            for i in range(len(nodes)):
                
                k = i + 1
                for j in range(k,len(nodes)):
            
                    if nodes[i] !=nodes[j]:
                        
                        temp = distance[nodes[i],nodes[j]] - average_cluster[clust]
                        sum_var  = sum_var + pow(temp,2) 
                        
                        count = count + 1 
                        #print(nodes[i],nodes[j])
                    
        
            variance_cluster[clust] = sum_var/count
        else:
            average_cluster[clust] = 1
            variance_cluster[clust] = 0
            
            
    return average_cluster,variance_cluster

def Transform_labels_2dList_to_intArray(distance,lab):
    labels = np.zeros(len(distance),dtype="int")
    n_clusters = len(lab[:])
    for k in range(0,n_clusters):
        
        for l in range(0,len(lab[k][:])) :
            
            # k is the cluster
            # l is the vertex node 
            i_index = lab[k][l]
            labels[i_index] = k
            
    return labels

def Plot_graph_visual_style(clust,name_of_graph):   
    
    visual_style = dict()
    visual_style["bbox"] = (500,500)
    #visual_style["vertex_label"] = labels
    ig.plot(clust,name_of_graph,mark_groups = False,**visual_style)
    return

def add_weights_as_list(g,distance):
    
    '''
        Input : graph 
                distance/similarity matrix
        -----------------------------
        Output: graph with attribute weights completed
    
    '''
    

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

def Plot_Graph(graph,label,name_of_graph,layout):
    #graph.vs["label"] = label
    color_dict = []
    values, counts = np.unique(label, return_counts=True)
    n = len(values)+1 #clusters start from 0 or 1!
    for i in range(n):
        color_dict.append('#%06X' % randint(0, 0xFFFFFF))
    
    
    graph.vs["color"] = [color_dict[l] for l in label]
    ig.plot(graph, name_of_graph, layout=layout, bbox=(500, 500), margin=50, inline='None')

    return 

def Plot_Distribution(distance):
    '''
        Input : similarity matrix
        -------------------------
        Output : figure object
    '''
    A = distance.flatten()
    
    A = np.log10(A*10)
    # t = sns.distplot(A,hist=True, kde=True,
    #                  bins=int(300 / 5), color="darkblue",
    #                  hist_kws={'edgecolor': 'black'},
    #                  kde_kws={'linewidth': 2})
    t = sns.histplot(x=A, kde=True,stat="probability",
                     bins=int(300 / 5))

    figure = t.get_figure()
    return figure,A

def calculate_distance(data):
    
    '''
        Input  : a matrix with string sequences 
        --------------------------------------
        Output : a similarity matrix normalized to [0,1]
    '''
    distance = np.zeros((len(data), len(data)))

    
    for i in range(0, len(distance)):
        start = i 

        for j in range(start,len(distance)):
            string1 = data.iloc[i, 1]  # take the second column
            string2 = data.iloc[j, 1]
            temp_distance = 1 - editdistance.eval(string1, string2) / len(string1)
            # fast implementation of Levenshtein
            # https://github.com/roy-ht/editdistance
            temp_distance = format(temp_distance, '.3f')
            distance[i, j] = temp_distance
            distance[j, i] = temp_distance

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

    graph_threshold = ig.Graph.Adjacency((distance > threshold).tolist())  
    # create a graph from adjacency matrix

    return graph_threshold


def normalized_weighted_graph(distance):
    
    ''' 
        Input : similarity matrix
        -------------------------
        Output: normalized_vertex_graph
    '''
    size = distance.shape[0]
    #global adj_matrix_boolean
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
    
    SAMPLES = ["B1",'S1',"T1"]
    
    
    for SAMPLE in SAMPLES:
        FILE = '2_IMGT-gapped-nt-sequences_' + SAMPLE + '.txt'
    
    
        COLUMN = 'V-D-J-REGION'
        SEQUENCE_NUMBER = 'Sequence number'
        
        path = os.getcwd()
        e = datetime.datetime.now()
        Hour_time_name = str(e.hour) + "_" + str(e.minute) + "_" +  str(e.second)
        
        Day_month = str(e.day) + "_"  + str(e.month) + "_"
        
    
        OUTPUT_FOLDER = SAMPLE
        OUTPUT_FOLDER_IMAGES = SAMPLE + "/IMAGES/"
    
        #os.makedirs(OUTPUT_FOLDER)
        f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "w")
        
        
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        data = pd.read_csv(os.path.join(__location__, FILE), sep='\t')
        filtered_data = data[[SEQUENCE_NUMBER, COLUMN]].sample(n=500).dropna()# drop nan # take 500 random samples
        
    
        ''' this block calculates distance '''
        start = time.time()
        distance = calculate_distance(filtered_data)
        
        
        ''' ***************************    '''
        
        
        ''' this block calculates the distribution of similarities '''
        
    
        figure,A = Plot_Distribution(distance)
        
        figure.savefig(OUTPUT_FOLDER +'/distribution_of_similarities.jpg', dpi=400)
        
        ''' ***************************    '''
    
        #--------------------------------------------------------
        layout = "kk"
        
        ''' This block produces the normalized graph '''
    
        graph_normalized = normalized_weighted_graph(distance)
            
        ig.plot(graph_normalized,OUTPUT_FOLDER + "graph_normalized.png", 
                                    layout=layout, bbox=(500, 500), margin=50, inline='None')
        
        ''' ***************************    '''
    
        
        temp_graph          = graph_normalized # check if take normalized or graph threshold
        adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
        threshold_neighbors = 300
        graph_shared,label  = shared_nearest_neighbors(adj_matrix_boolean,threshold_neighbors)
            
            
        Plot_Graph(graph_shared,label,OUTPUT_FOLDER + "graph_normalized_shared.png",layout)
       
    
        ''' this block produces the analysis for each threshold '''
        
        threshold = np.array([0.70,0.75,0.80,0.85,0.90,0.95])
        
        for i in range(0,len(threshold[:])):
            
            
            graph_threshold = Choose_Threshold_Graph(distance,  threshold[i])
    
            threshold_name  = "threshold_" + str(threshold[i]) + "_"
        
            # ig.plot(graph_threshold, OUTPUT_FOLDER + threshold_name + "graph_threshold.png",
            #                         layout=layout, bbox=(500, 500), margin=50, inline='None')
        
            #--------------------------------------------------------
            temp_graph          = graph_threshold # check if take normalized or graph threshold
            adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
            threshold_neighbors = 300
            graph_shared,label  = shared_nearest_neighbors(adj_matrix_boolean,threshold_neighbors)
            # Plot_Graph(graph_shared,label,OUTPUT_FOLDER + threshold_name + "graph_shared.png",layout)
    
            #--------------------------------------------------------
            
            
            # g = graph_threshold
            # g = add_weights_as_list(g,distance)
            # wtrap = g.community_walktrap(weights = g.es['weight'],steps = 30)
            # lab =np.array(wtrap.as_clustering()).tolist()
            # clust_wtrap =  wtrap.as_clustering()
            # labels_wtrap = Transform_labels_2dList_to_intArray(distance,lab)
            # avg_clust_wtrap,var_clust_wtrap =Calculate_distance_metric_inside_cluster(distance,clust_wtrap)
    
            # Plot_Graph(g,labels_wtrap,OUTPUT_FOLDER + threshold_name + "clustering_wtrap_layout_kk.png",layout)
            # Plot_graph_visual_style(clust_wtrap,OUTPUT_FOLDER + threshold_name + "clustering_wtrap_visual_style.png")
    
            
            g = graph_threshold
            g.to_undirected()
            g = add_weights_as_list(g,distance)
            fastgreedy = g.community_fastgreedy(weights = g.es['weight'])    
            lab =np.array(fastgreedy.as_clustering()).tolist()
            clust_fastgreedy =  fastgreedy.as_clustering()
            labels_fastgreedy = Transform_labels_2dList_to_intArray(distance,lab)
            avg_clust_fastgreedy,var_clust_fastgreedy = Calculate_distance_metric_inside_cluster(distance, clust_fastgreedy)
            
            Plot_Graph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_layout_kk.png",layout)
            Plot_graph_visual_style(clust_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_visual_style.png")
        
        
            # layout = g.layout_drl(weights=g.es["weight"], fixed=None, seed=None, options=None, dim=2)
        
            # Plot_Graph(g,labels_fastgreedy,OUTPUT_FOLDER + threshold_name + "clustering_test.png",layout)
        
            f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "a")
            f_m.write("modularity for fast_greedy " +  threshold_name + "      " + str(clust_fastgreedy.modularity) +" \n")
            #f_m.write("average clust  similarity " +  threshold_name + str(np.round(avg_clust_fastgreedy,decimals=4)) +" \n")
            #f_m.write("variance clust similarity  " +  threshold_name + str(np.round(var_clust_fastgreedy,decimals=4))+" \n")
            f_m.write("--------------------------------------------------------------- \n \n")
        
        
        
        end = time.time()
        elapsed_time = format(end - start, '.3f')
        print('The execution time is {0} seconds\n'.format(elapsed_time))
        
        
        A = distance.flatten()
        mean_sim  = statistics.fmean(A)
        median_sim = statistics.median(A)
        pstdev_sim = statistics.pstdev(data = A , mu = mean_sim)
        
        percentile_50_sim = np.percentile(A,50)   
        percentile_75_sim = np.percentile(A,75)
        
        f = open(path + "\\" + OUTPUT_FOLDER + "\\similarity_metrics_" + SAMPLE + ".txt", "w")
        
        f.write("Some similarity metrics here \n---------------------------- \n")
        f.write("mean_sim  :  " + str(mean_sim) + " \n")
        f.write("median_sim  :  " + str(median_sim) +" \n")
        f.write("pstdev_sim  :  " + str(pstdev_sim) + "\n")
        f.write("percentile-50_sim  :  " + str(percentile_50_sim) + "\n")
        f.write("percentile-75_sim  :  " + str(percentile_75_sim) + "\n")
    
    
        f_m.close()
        f.close()
        
        
         
        # g = graph_threshold
        # g.to_undirected()
        # g = add_weights_as_list(g,distance)
        # clust_lead_eig = g.community_leading_eigenvector(clusters = 5 ,weights = g.es['weight'])    
    
    
        # g = graph_threshold
        # g.to_undirected()
        # g = add_weights_as_list(g,distance)
        # clust_optimal= g.community_optimal_modularity(weights = g.es['weight'])    
    
        # g = graph_threshold
        # g.to_undirected()
        # g = add_weights_as_list(g,distance)
        # clust_multi_level= g.community_multilevel(weights = g.es['weight'],return_levels=True) 
        #10 clusters
    
    
        #matrix_to_spectral = graph_normalized.get_adjacency()[:, :]  # int 0 1
        # matrix_to_spectral = graph_threshold.get_adjacency()[:, :]  # int 0 1
        
        # matrix_to_spectral = np.array(matrix_to_spectral.data)  # int
        
        # d = np.diag(matrix_to_spectral.sum(axis=1))
        
        # l = np.array(d - matrix_to_spectral)
        # eig_L, eig_vec_L = np.linalg.eig(l)
        
        # i = np.where(eig_L < 0.1)
        
        #sc               = SpectralClustering(5, affinity='precomputed', n_init=100,assign_labels='discretize')   
        
        #label           = sc.fit_predict(adjacency_matrix)
                        
        end = time.time()    

    