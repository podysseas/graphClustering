



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
    b = filtered_data['labels']
    prep_heatmap = {gene: np.array(a), 'labels':b}
    
    df = pd.DataFrame(prep_heatmap)
    
    clusters_freq = df.groupby(['labels',gene]).size().reset_index(name="freq")
    pivoted_clusters = clusters_freq.pivot_table(index=gene, 
              columns='labels', 
              values='freq',
              fill_value=0).unstack().to_frame()
    pivoted_clusters.reset_index(inplace=True)  
    
    pivoted_clusters.columns = ['labels', gene, 'count']
    pivoted_clusters = pivoted_clusters.groupby(['labels',gene]).agg({'count': 'sum'})
    pivoted_clusters = pivoted_clusters.groupby(level=0).apply(lambda x:
                                                 100 * x / float(x.sum()))
    pivoted_clusters.reset_index(inplace=True)  
    
    for item in all_genes:
        
        total_df[item] = np.around(np.array(pivoted_clusters[pivoted_clusters[gene]==item]['count'],dtype=float),decimals=2)
    
    t = sns.heatmap(total_df.T,vmax=100 )
    for item in t.get_yticklabels():
        item.set_rotation(45)
        
    figure = t.get_figure()
    
    figure.savefig(pathToSave)
    plt.close(figure)
    #pivoted_clusters['count']  = np.array(pivoted_clusters['count'],dtype=int)
    
    return