import numpy as np
import statistics
import igraph as ig
import editdistance

def TransformLabels2dListToIntArray(distance,lab):
    labels = np.zeros(len(distance),dtype="int")
    n_clusters = len(lab[:])
    for k in range(0,n_clusters):
        
        for l in range(0,len(lab[k][:])) :
            
            # k is the cluster
            # l is the vertex node 
            i_index = lab[k][l]
            labels[i_index] = k
            
    return labels

def normalizedWeightedGraph(distance):
    
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

def addWeightsAsList(g,distance):
    
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



def calculateDistance(data):
    
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
            if i == j :
                temp_distance = 0
            else:
                temp_distance = 1 - editdistance.eval(string1, string2) / len(string1)
            # fast implementation of Levenshtein
            # https://github.com/roy-ht/editdistance
            temp_distance = format(temp_distance, '.3f')
            distance[i, j] = temp_distance
            distance[j, i] = temp_distance

    return distance


def ChooseThresholdGraph(distance, threshold):
    
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

def cutEdgesDynamicThreshold(distance):
    
    size = len(distance)
    adj_matrix_boolean = np.ones((size, size), dtype=int)

    
    for col in range(0,len(distance)):
        
        distance[:,col] = np.floor(distance[:,col]*100)/100 # floor to second decimal
        
        
    max_of_columns = np.amax(distance, axis=1)   # Maxima along the second axis | columns
    
    for i in range(0, size):
        
        temp_threshold = max_of_columns[i]
        
        # vertex_threshold[i] = format(temp_threshold, '.4f')

        adj_matrix_boolean[i, :] = distance[i, :] >= temp_threshold

    graph_dynamic = ig.Graph.Adjacency(adj_matrix_boolean.tolist())  # create a graph from adjacency matrix

    return graph_dynamic