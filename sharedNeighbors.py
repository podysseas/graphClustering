
import numpy as np

import igraph as ig

#im
def sharedNearestNeighbors(adj_matrix_boolean,threshold_neighbors):
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
    
    
    return graph_shared,label

