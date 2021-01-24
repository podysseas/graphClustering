"""
Created on Sat Sep 19 17:10:24 2020

@author: ody
"""

"""

"""
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import os
import pandas as pd
import time
import igraph as ig


#import matplotlib.pyplot as plt
import markov_clustering as mc

from datasetMetrics import calculateMetrics,   CalculateDistanceMetricInsideCluster
from helpFunctions import ChooseThresholdGraph,  calculateDistance,normalizedWeightedGraph, TransformLabels2dListToIntArray, addWeightsAsList
from plotFunctions import PlotDistribution ,  PlotGraph,  PlotGraphVisualStyle
from sharedNeighbors import  sharedNearestNeighbors
from heatmap import heatmap
if __name__ == "__main__":
    
    SAMPLES = ["T1","B1","S1"]
    
    gene =  "J-GENE and allele" 
    
    
    start = time.time()

    for SAMPLE in SAMPLES:
        
        FILE = '2_IMGT-gapped-nt-sequences_' + SAMPLE + '.txt'
        COLUMN = 'V-D-J-REGION'
        SEQUENCE_NUMBER = 'Sequence number'
        
        path = os.getcwd()
        
        print("#----------------------------------#\n  SAMPLE : =  " + str(SAMPLE) + "\n\n")
        OUTPUT_FOLDER = SAMPLE
        OUTPUT_FOLDER_IMAGES = SAMPLE + "/IMAGES/"
    
        
        f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "w")
        f_m.write(" ----------  Modularities  ------------ \n\n")
        
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
        data = pd.read_csv(os.path.join(__location__, FILE), sep='\t')
        data.columns = ['Sequence number', 'Sequence ID', 'V-DOMAIN Functionality',
                           'V-GENE and allele', 'J-GENE and allele', 'D-GENE and allele',
                           'V-D-J-REGION', 'V-J-REGION', 'V-REGION', 'FR1-IMGT', 'CDR1-IMGT',
                           'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT', 'JUNCTION',
                           'J-REGION', 'FR4-IMGT', 'Unnamed: 18']
        filtered_data = data[[SEQUENCE_NUMBER, COLUMN]].sample(n=500 , random_state=42).dropna()# drop nan # take 500 random samples
        
        filtered_data = data[[SEQUENCE_NUMBER, COLUMN,gene]].sample(n=500).dropna()# drop nan # take 500 random samples
        all_genes = filtered_data[gene].unique()

   
        ''' this block calculates distance '''
        distance = calculateDistance(filtered_data)
        
        
        
        
        ''' this block calculates some metrics for this datasets '''
        calculateMetrics(distance,path,OUTPUT_FOLDER,SAMPLE)        
        
        
        ''' this block calculates the distribution of similarities '''

        figure,A = PlotDistribution(distance)
        
        figure.savefig(OUTPUT_FOLDER +'/distribution_of_similarities.jpg', dpi=400)
    
        #--------------------------------------------------------
        layout = "kk"
        
        ''' This block produces the normalized graph '''
    
        graph_normalized = normalizedWeightedGraph(distance)
            
        ig.plot(graph_normalized,OUTPUT_FOLDER + "graph_normalized.png", 
                                    layout=layout, bbox=(500, 500), margin=50, inline='None')
        
    
        
        # temp_graph          = graph_normalized # check if take normalized or graph threshold
        # adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
        # threshold_neighbors = 300
        # graph_shared,label  = shared_nearest_neighbors(adj_matrix_boolean,threshold_neighbors)
            
            
        # PlotGraph(graph_shared,label,OUTPUT_FOLDER + "graph_normalized_shared.png",layout)
       
    
        ''' this block produces the analysis for each threshold '''
        
        threshold = np.array([0.70,0.80,0.85,0.9])
        
        for i in range(0,len(threshold)):
            
            
            graph_threshold = ChooseThresholdGraph(distance,  threshold[i])
    
            threshold_name  = "threshold_" + str(threshold[i]) + "_"
        
            # ig.plot(graph_threshold, OUTPUT_FOLDER + threshold_name + "graph_threshold.png",
            #                         layout=layout, bbox=(500, 500), margin=50, inline='None')
        
            #--------------------------------------------------------
            temp_graph          = graph_threshold # check if take normalized or graph threshold
            adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
            threshold_neighbors = 300
            graph_shared,label  = sharedNearestNeighbors(adj_matrix_boolean,threshold_neighbors)
            # PlotGraph(graph_shared,label,OUTPUT_FOLDER + threshold_name + "graph_shared.png",layout)
    
            #--------------------------------------------------------
            
            
            # g = graph_threshold
            # g = add_weights_as_list(g,distance)
            # wtrap = g.community_walktrap(weights = g.es['weight'],steps = 30)
            # lab =np.array(wtrap.as_clustering()).tolist()
            # clust_wtrap =  wtrap.as_clustering()
            # labels_wtrap = TransformLabels2dListToIntArray(distance,lab)
            # avg_clust_wtrap,var_clust_wtrap =Calculate_distance_metric_inside_cluster(distance,clust_wtrap)
    
            # Plot_Graph(g,labels_wtrap,OUTPUT_FOLDER + threshold_name + "clustering_wtrap_layout_kk.png",layout)
            # Plot_graph_visual_style(clust_wtrap,OUTPUT_FOLDER + threshold_name + "clustering_wtrap_visual_style.png")
    
            
            g = graph_threshold
            g.to_undirected()
            g = addWeightsAsList(g,distance)
            fastgreedy = g.community_fastgreedy(weights = g.es['weight'])    
            lab =np.array(fastgreedy.as_clustering()).tolist()
            clust_fastgreedy =  fastgreedy.as_clustering()
            labels_fastgreedy = TransformLabels2dListToIntArray(distance,lab)
            avg_clust_fastgreedy,var_clust_fastgreedy = CalculateDistanceMetricInsideCluster(distance, clust_fastgreedy)
            
            PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_layout_kk.png",layout)
            PlotGraphVisualStyle(clust_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_visual_style.png")
        
        
            # layout = g.layout_drl(weights=g.es["weight"], fixed=None, seed=None, options=None, dim=2)
        
            # Plot_Graph(g,labels_fastgreedy,OUTPUT_FOLDER + threshold_name + "clustering_test.png",layout)
        
            
            
            ''' ########################## '''
            pathToSave = OUTPUT_FOLDER +'/'+ "heatmap_" + threshold_name + '.jpg'
            heatmap(gene,filtered_data, labels_fastgreedy ,all_genes, pathToSave)            
            ''' '''
            
            g = graph_threshold
            A = np.array(graph_threshold.get_adjacency().data)
            g.to_undirected()

            
            inflation = 2
            result = mc.run_mcl(A, inflation=inflation)
            clusters = mc.get_clusters(result)
            
            final_clusters  = np.zeros(len(labels_fastgreedy))
            for j in range(0,len(clusters)):
                np.add.at(final_clusters,np.array(clusters[j]),np.full(len(clusters[j]),j))
       
            #Q = mc.modularity(matrix=result, clusters=clusters)
            g.vs['group']= final_clusters
            Q = g.modularity(g.vs["group"])
    
            
    
            f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "a")

            f_m.write("Threshold :  " + str(threshold[i]) + "\n\n")
            
            f_m.write("MCL algorithm - modularity:   " + str(Q) + "\n\n")
            
            f_m.write("Fast Greedy algorithm:   " + str(clust_fastgreedy.modularity) +" \n\n\n")
            
            
            #f_m.write("average clust  similarity " +  threshold_name + str(np.round(avg_clust_fastgreedy,decimals=4)) +" \n")
            #f_m.write("variance clust similarity  " +  threshold_name + str(np.round(var_clust_fastgreedy,decimals=4))+" \n")
            f_m.write("--------------------------------------------------------------- \n \n")
        



        g = graph_normalized
        g.to_undirected()
        A = np.array(graph_normalized.get_adjacency().data)
        inflation = 2
        result = mc.run_mcl(A, inflation=inflation)
        clusters = mc.get_clusters(result)
        
        final_clusters  = np.zeros(len(labels_fastgreedy))
        for j in range(0,len(clusters)):
            np.add.at(final_clusters,np.array(clusters[j]),np.full(len(clusters[j]),j))
   
        #Q = mc.modularity(matrix=result, clusters=clusters)
        g.vs['group']= final_clusters
        Q = g.modularity(g.vs["group"])
        
        f_m.write("Graph normalized  : ")
        f_m.write("MCL algorithm - modularity:   " + str(Q) + "\n\n")

        f_m.close()

        end = time.time()
        elapsed_time = format(end - start, '.3f')
        print('The execution time is {0} seconds\n'.format(elapsed_time))









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
                        

    