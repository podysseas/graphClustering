
import numpy as np
import seaborn as sns
import igraph as ig

from random import randint



def PlotGraph(graph,label,name_of_graph,layout):
    #graph.vs["label"] = label
    color_dict = []
    values, counts = np.unique(label, return_counts=True)
    n = len(values)+1 #clusters start from 0 or 1!
    for i in range(n):
        color_dict.append('#%06X' % randint(0, 0xFFFFFF))
    
    
    graph.vs["color"] = [color_dict[l] for l in label]
    ig.plot(graph, name_of_graph, layout=layout, bbox=(500, 500), margin=50, inline='None')

    return 

def PlotDistribution(distance):
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
    t.set_xlim(0.65,1)
    figure = t.get_figure()
    return figure,A

def PlotGraphVisualStyle(clust,name_of_graph):   
    
    visual_style = dict()
    visual_style["bbox"] = (500,500)
    #visual_style["vertex_label"] = labels
    ig.plot(clust,name_of_graph,mark_groups = False,**visual_style)
    return


