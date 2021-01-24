



import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
def heatmap(gene,filtered_data,labels,all_genes,pathToSave):
    
    filtered_data['labels'] = labels
    
    total_df = pd.DataFrame(columns=all_genes)
    total_df.reset_index( inplace=False)
    
       
    
    '''This block produces the HEATMAP'''
    a = filtered_data[gene]
    b = labels[labels]
    prep_heatmap = {'J-GENE and allele': np.array(a), 'labels':b}
    
    df = pd.DataFrame(prep_heatmap)
    
    clusters_freq = df.groupby(['labels',gene]).size().reset_index(name="freq")
    
    pivoted_clusters = clusters_freq.pivot_table(index=gene, 
              columns='labels', 
              values='freq',
              fill_value=0).unstack().to_frame()
    pivoted_clusters.reset_index(inplace=True)  
    pivoted_clusters.columns = ['labels', 'J-GENE and allele', 'count']
    
    for item in all_genes:
        
        total_df[item] = np.array(pivoted_clusters[pivoted_clusters[gene]==item]['count'],dtype=int)
    
    t = sns.heatmap(total_df.T)
    figure = t.get_figure()
    
    figure.savefig(pathToSave)
    plt.close(figure)

    return