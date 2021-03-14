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

#from test import cutEdgesDynamicThreshold
from natsort import natsorted
from datasetMetrics import calculateMetrics,   CalculateDistanceMetricInsideCluster
from helpFunctions import ChooseThresholdGraph,  calculateDistance,normalizedWeightedGraph,TransformLabels2dListToIntArray, addWeightsAsList,cutEdgesDynamicThreshold
from plotFunctions import PlotDistribution ,  PlotGraph #,  PlotGraphVisualStyle
from sharedNeighbors import  sharedNearestNeighbors
from heatmap import heatmap

if __name__ == "__main__":
    
    SAMPLES = ["T1"]
    
    genes =  ['J-GENE and allele','V-GENE and allele','D-GENE and allele']
    
    
    start = time.time()

    for SAMPLE in SAMPLES:
        
        
            FILE = '2_IMGT-gapped-nt-sequences_' + SAMPLE + '.txt'
            COLUMN = 'V-D-J-REGION'
            SEQUENCE_NUMBER = 'Sequence number'
            
            path = os.getcwd()
            
            print("#----------------------------------#\n  SAMPLE : =  " + str(SAMPLE) + "\n\n")
            OUTPUT_FOLDER = SAMPLE
            OUTPUT_FOLDER_IMAGES = SAMPLE + "/IMAGES/"
            OUTPUT_FOLDER_GENES = SAMPLE + "/GENES/"
            OUTPUT_FOLDER_HEATMAP = SAMPLE + "/HEATMAP/"

            
            f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "w")
            f_m.write(" ----------  Modularities  ------------ \n\n")
            
            __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
            data = pd.read_csv(os.path.join(__location__, FILE), sep='\t')
            data.columns = ['Sequence number', 'Sequence ID', 'V-DOMAIN Functionality',
                               'V-GENE and allele', 'J-GENE and allele', 'D-GENE and allele',
                               'V-D-J-REGION', 'V-J-REGION', 'V-REGION', 'FR1-IMGT', 'CDR1-IMGT',
                               'FR2-IMGT', 'CDR2-IMGT', 'FR3-IMGT', 'CDR3-IMGT', 'JUNCTION',
                               'J-REGION', 'FR4-IMGT', 'Unnamed: 18']
            #filtered_data = data[[SEQUENCE_NUMBER, COLUMN]].sample(n=500 , random_state=42).dropna()# drop nan # take 500 random samples
            
            filtered_data = data[[SEQUENCE_NUMBER, COLUMN,genes[0],genes[1],genes[2]]].sample(n=500,random_state=50).dropna()# drop nan # take 500 random samples

            # preprocess the column for the gene 
           
               
            
    
       
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
            
            graph_dynamic   = cutEdgesDynamicThreshold(distance)
        
            
            # temp_graph          = graph_normalized # check if take normalized or graph threshold
            # adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
            # threshold_neighbors = 300
            # graph_shared,label  = shared_nearest_neighbors(adj_matrix_boolean,threshold_neighbors)
                
                
            # PlotGraph(graph_shared,label,OUTPUT_FOLDER + "graph_normalized_shared.png",layout)
           
        
            ''' this block produces the analysis for each threshold '''
            
            threshold = np.array([0.70,0.80,0.85,0.9])
            
            for th in range(0,len(threshold)):
                
                
                graph_threshold = ChooseThresholdGraph(distance,  threshold[th])
        
                threshold_name  = "threshold_" + str(threshold[th]) + "_"
            
                # ig.plot(graph_threshold, OUTPUT_FOLDER + threshold_name + "graph_threshold.png",
                #                         layout=layout, bbox=(500, 500), margin=50, inline='None')
            
                #--------------------------------------------------------
                # temp_graph          = graph_threshold # check if take normalized or graph threshold
                # adj_matrix_boolean  = np.array(temp_graph.get_adjacency().data)
                # threshold_neighbors = 300
                # graph_shared,label  = sharedNearestNeighbors(adj_matrix_boolean,threshold_neighbors)
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
        
                # graph_threshold.to_undirected()
                # graph_threshold = addWeightsAsList(graph_threshold,distance)
                # leiden = graph_threshold.community_leiden(objective_function="modularity", 
                #                                             weights=graph_threshold.es['weight'], 
                #                                             resolution_parameter=0.4, beta=0.5, 
                #                                             initial_membership=None, n_iterations=3, node_weights=None)
                
                # g = graph_threshold
                graph_threshold.to_undirected()
                graph_threshold = addWeightsAsList(graph_threshold,distance)
                fastgreedy = graph_threshold.community_fastgreedy(weights = graph_threshold.es['weight'])    
                lab =np.array(fastgreedy.as_clustering()).tolist()
                clust_fastgreedy =  fastgreedy.as_clustering()
                labels_fastgreedy = TransformLabels2dListToIntArray(distance,lab)
                avg_clust_fastgreedy,var_clust_fastgreedy = CalculateDistanceMetricInsideCluster(distance, clust_fastgreedy)
                # fruchterman
                # PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_layout_kk.png","kk")
                
                
                # PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_layout_fruchterman.png","fr")
    
                # PlotGraphVisualStyle(clust_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_visual_style.png")
            
                # PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_layout_lgl.png","lgl")
    
                # layout = g.layout_drl(weights=g.es["weight"], fixed=None, seed=None, options=None, dim=2)
            
                # PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES + threshold_name + "clustering_fastgreedy_layout_kk.png","kk")
                
                # PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES +  threshold_name + "clustering_fastgreedy_layout_rt_circular.png","rt_circular")
    
                # the layout_mds 
                
                mds_layout = graph_threshold.layout_mds(dist=distance,dim=2)
                PlotGraph(graph_threshold,labels_fastgreedy,OUTPUT_FOLDER_IMAGES + threshold_name + "clustering_fastgreedy_layout_mds.png",mds_layout)
    
                # PlotGraph(g,labels_fastgreedy,OUTPUT_FOLDER_IMAGES + threshold_name + "clustering_fastgreedy_layout_auto.png","auto")
    
                
    
                ''' ########################## '''
                for gene in genes:
                    
                    genes_data = filtered_data[[SEQUENCE_NUMBER, COLUMN,gene]]# drop nan # take 500 random samples

                    

                    temp_gene = np.array(genes_data[gene])
                    if SAMPLE=="B1":
                        temp_gene[159] = "Homsap IGHJ3"
                        temp_gene[396] = "Homsap IGHJ4"
            
                    for ge in range(len(temp_gene)):
                        
                        temp_gene[ge]  = str(temp_gene[ge]).split("*")[0].split(" ")[1]
                         
                    genes_data[gene] = temp_gene
                    
                    all_genes = natsorted( genes_data[gene].unique() )
                    
                    pathToSave = OUTPUT_FOLDER_HEATMAP +'/'+ "heatmap_" + threshold_name + "_" + str(gene.split("-")[0]) + '_.jpg'
                    heatmap(gene,genes_data, labels_fastgreedy ,all_genes, pathToSave) 
               
                
                
                    # step 1 : map all genes to one number
                    # step 2 : produce the labels array according to the dictionary d
                    # step 3 : make it a function 
                    #        - input  : genes_column , all_genes ( not necessary - maybe with unique )
                    #        - output : returns the array with the labels_genes
                    # step 4 : plot the graph with these labels
                    
                    # ------ plot the same graph but now labels are from genes ------
                    d = dict([(y,x+1) for x,y in enumerate(sorted(set(all_genes)))])
        
                    
                    labels_genes = np.array([],dtype=int)
                    for my_g in genes_data[gene]:
                        
                        # dict['Age']:  8
                        
                        labels_genes = np.append( labels_genes , d[my_g] )
                    
                    # pathToSave = OUTPUT_FOLDER +'/'+ "label_genes_" + threshold_name + 'fruchterman.png'
             
                    # PlotGraph(g,labels_genes,pathToSave,"fr")
        
                    pathToSave = OUTPUT_FOLDER_GENES +'/'+ "label_genes_" + threshold_name + str(gene) + '_mds.png'
             
                    PlotGraph(graph_threshold,labels_genes,pathToSave,mds_layout)
                
                
                # --------------------------
                # try the mcl algorithm 
                
                
                # g = graph_threshold
                A = np.array(graph_threshold.get_adjacency().data)
                graph_threshold.to_undirected()
    
                
                inflation = 2
                result = mc.run_mcl(A, inflation=inflation)
                clusters = mc.get_clusters(result)
                
                final_clusters  = np.zeros(len(labels_fastgreedy))
                for j in range(0,len(clusters)):
                    np.add.at(final_clusters,np.array(clusters[j]),np.full(len(clusters[j]),j))
           
                #Q = mc.modularity(matrix=result, clusters=clusters)
                graph_threshold.vs['group']= final_clusters
                Q = graph_threshold.modularity(graph_threshold.vs["group"])
        
                # -------------------------------
        
                f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "a")
    
                f_m.write("Threshold :  " + str(threshold[th]) + "\n\n")
                
                f_m.write("MCL algorithm - modularity:   " + str(Q) + "\n\n")
                
                f_m.write("Fast Greedy algorithm:   " + str(clust_fastgreedy.modularity) +" \n\n\n")
                
                
                #f_m.write("average clust  similarity " +  threshold_name + str(np.round(avg_clust_fastgreedy,decimals=4)) +" \n")
                #f_m.write("variance clust similarity  " +  threshold_name + str(np.round(var_clust_fastgreedy,decimals=4))+" \n")
                f_m.write("--------------------------------------------------------------- \n \n")
        
        
            # g = graph_normalized
            # g.to_undirected()
            # A = np.array(graph_normalized.get_adjacency().data)
            # inflation = 2
            # result = mc.run_mcl(A, inflation=inflation)
            # clusters = mc.get_clusters(result)
        
        
        
            # final_clusters  = np.zeros(len(labels_fastgreedy),dtype=int)
            # for j in range(0,len(clusters)):
            #     np.add.at(final_clusters,np.array(clusters[j]),np.full(len(clusters[j]),j))
       
            # #Q = mc.modularity(matrix=result, clusters=clusters)
            # g.vs['group']= final_clusters
            # Q = g.modularity(g.vs["group"])
            
            # f_m.write("Graph normalized  : ")
            # f_m.write("MCL algorithm - modularity:   " + str(Q) + "\n\n")
            
            graph_normalized.to_undirected()
            graph_normalized = addWeightsAsList(graph_normalized,distance)
            fastgreedy = graph_normalized.community_fastgreedy(weights = graph_normalized.es['weight'])    
            lab =np.array(fastgreedy.as_clustering()).tolist()
            clust_fastgreedy =  fastgreedy.as_clustering()
            labels_fastgreedy = TransformLabels2dListToIntArray(distance,lab)
            avg_clust_fastgreedy,var_clust_fastgreedy = CalculateDistanceMetricInsideCluster(distance, clust_fastgreedy)  
            
            mds_layout = graph_normalized.layout_mds(dist=distance,dim=2)
            PlotGraph(graph_normalized,labels_fastgreedy,OUTPUT_FOLDER_IMAGES + "graph_normalized_" + "clustering_fastgreedy_layout_mds.png",mds_layout)
    
            
            for gene in genes:
                        
                        genes_data = filtered_data[[SEQUENCE_NUMBER, COLUMN,gene]]# drop nan # take 500 random samples

                       
    
                        temp_gene = np.array(genes_data[gene])
                        if SAMPLE=="B1":
                            temp_gene[159] = "Homsap IGHJ3"
                            temp_gene[396] = "Homsap IGHJ4"
                
                        for ge in range(len(temp_gene)):
                            
                            temp_gene[ge]  = str(temp_gene[ge]).split("*")[0].split(" ")[1]
                             
                        genes_data[gene] = temp_gene
                        all_genes = natsorted( genes_data[gene].unique() )
                         
                        pathToSave = OUTPUT_FOLDER_HEATMAP +'/'+ "heatmap_" + "graph_normalized_" + str(gene.split("-")[0]) + '_.jpg'
                        heatmap(gene,genes_data, labels_fastgreedy ,all_genes, pathToSave) 
                        
                        d = dict([(y,x+1) for x,y in enumerate(sorted(set(all_genes)))])
            
                        
                        labels_genes = np.array([],dtype=int)
                        for my_g in genes_data[gene]:
                            
                            # dict['Age']:  8
                            
                            labels_genes = np.append( labels_genes , d[my_g] )
                        
                        # pathToSave = OUTPUT_FOLDER +'/'+ "label_genes_" + threshold_name + 'fruchterman.png'
                 
                        # PlotGraph(g,labels_genes,pathToSave,"fr")
        
                        pathToSave = OUTPUT_FOLDER_GENES +'/'+ "label_genes_" + "graph_normalized_" + str(gene) + '_mds.png'
                 
                        PlotGraph(graph_normalized,labels_genes,pathToSave,mds_layout)
            
            # here g is graph normalized
            # g = graph_normalized
            
                    
            # -----------------------------------------------
            # --------------------------------------------------        
            # g = graph_dynamic
            # graph_dynamic.to_undirected()
            # A = np.array(graph_dynamic.get_adjacency().data)
            # inflation = 2
            # result = mc.run_mcl(A, inflation=inflation)
            # clusters = mc.get_clusters(result)
            
            # final_clusters  = np.zeros(len(labels_fastgreedy))
            # for j in range(0,len(clusters)):
            #     np.add.at(final_clusters,np.array(clusters[j]),np.full(len(clusters[j]),j))
       
            # #Q = mc.modularity(matrix=result, clusters=clusters)
            # graph_dynamic.vs['group']= final_clusters
            # Q = graph_dynamic.modularity(graph_dynamic.vs["group"])
            
            # f_m.write("Graph dynamic  : ")
            # f_m.write("MCL algorithm - modularity:   " + str(Q) + "\n\n")
    
            # here g is graph dynamic
            graph_dynamic.to_undirected()
            graph_dynamic = addWeightsAsList(graph_dynamic,distance)
            fastgreedy = graph_dynamic.community_fastgreedy(weights = graph_dynamic.es['weight'])    
            lab =np.array(fastgreedy.as_clustering()).tolist()
            clust_fastgreedy =  fastgreedy.as_clustering()
            labels_fastgreedy = TransformLabels2dListToIntArray(distance,lab)
            avg_clust_fastgreedy,var_clust_fastgreedy = CalculateDistanceMetricInsideCluster(distance, clust_fastgreedy)  
            
            mds_layout = graph_dynamic.layout_mds(dist=distance,dim=2)
            PlotGraph(graph_dynamic,labels_fastgreedy,OUTPUT_FOLDER_IMAGES + "graph_dynamic_" + "clustering_fastgreedy_layout_mds.png",mds_layout)
    
            for gene in genes:
                        
                        genes_data = filtered_data[[SEQUENCE_NUMBER, COLUMN,gene]]# drop nan # take 500 random samples

    
                        temp_gene = np.array(genes_data[gene])
                        if SAMPLE=="B1":
                            temp_gene[159] = "Homsap IGHJ3"
                            temp_gene[396] = "Homsap IGHJ4"
                
                        for ge in range(len(temp_gene)):
                            
                            temp_gene[ge]  = str(temp_gene[ge]).split("*")[0].split(" ")[1]
                             
                        genes_data[gene] = temp_gene
                        all_genes = natsorted( genes_data[gene].unique() )
                         
                        pathToSave = OUTPUT_FOLDER_HEATMAP +'/'+ "heatmap_" + "graph_dynamic_" + str(gene.split("-")[0]) + '_.jpg'
                        heatmap(gene,genes_data, labels_fastgreedy ,all_genes, pathToSave)
                        
                        d = dict([(y,x+1) for x,y in enumerate(sorted(set(all_genes)))])
            
                        
                        labels_genes = np.array([],dtype=int)
                        for my_g in genes_data[gene]:
                            
                            # dict['Age']:  8
                            
                            labels_genes = np.append( labels_genes , d[my_g] )
                        
                        # pathToSave = OUTPUT_FOLDER +'/'+ "label_genes_" + threshold_name + 'fruchterman.png'
                 
                        # PlotGraph(g,labels_genes,pathToSave,"fr")
        
                        pathToSave = OUTPUT_FOLDER_GENES +'/'+ "label_genes_" + "graph_dynamic_"+ str(gene) + '_mds.png'
                 
                        PlotGraph(graph_dynamic,labels_genes,pathToSave,mds_layout)
            f_m.close()

    end = time.time()
    elapsed_time = format(end - start, '.3f')
    print('The execution time is {0} seconds\n'.format(elapsed_time))























        # g = graph_normalized
        # g.to_undirected()
        # g = addWeightsAsList(g,distance)
        # clust_lead_eig = g.community_leading_eigenvector(clusters = 5 ,weights = g.es['weight'])    
    
    
        # g = graph_normalized
        # g.to_undirected()
        # g = addWeightsAsList(g,distance)
        # clust_multi_level= g.community_multilevel(weights = g.es['weight'],return_levels=True) 
        
        # f_m = open(path + "\\" + OUTPUT_FOLDER + "\\modularities_" + SAMPLE + ".txt", "a")
        # #f_m.write("Graph normalized  : ")
        # f_m.write("leading eigenvector  - modularity:   " + str(clust_lead_eig.modularity) + "\n\n")
        
        # #f_m.write("community optimal modularity   - modularity:   " + str(clust_optimal.modularity) + "\n\n")
        # f_m.write("multilevel  - modularity:   " + str(clust_multi_level[0].modularity) + "\n\n")

        # f_m.close()
    

                        

    
