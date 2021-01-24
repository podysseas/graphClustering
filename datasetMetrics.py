
import statistics
import numpy as np

def calculateMetrics(distance,path,OUTPUT_FOLDER,SAMPLE):
    
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
        f.close()
        return 
    
def CalculateDistanceMetricInsideCluster(distance,clusters):
    
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
