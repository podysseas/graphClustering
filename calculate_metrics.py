
import statistics
import numpy as np

def calculate_metrics(distance,path,OUTPUT_FOLDER,SAMPLE):
    
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