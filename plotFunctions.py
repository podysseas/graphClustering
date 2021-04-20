
import numpy as np
import seaborn as sns
import igraph as ig
import cairo
from random import randint
from igraph.drawing.text import TextDrawer


def PlotGraph(graph,label,name_of_graph,layout):
    #graph.vs["label"] = label
    color_dict = []
    values, counts = np.unique(label, return_counts=True)
    n = len(values)+1 #clusters start from 0 or 1!
    for i in range(n):
        color_dict.append('#%06X' % randint(0, 0xFFFFFF))
    
    
    graph.vs["color"] = [color_dict[l] for l in label]
    
    #bbox = BoundingBox(600, 600)

    myPlot = ig.plot(name_of_graph, bbox=(600,600), background="white")
    # myPlot = ig.plot(graph, name_of_graph, layout=layout, mark_groups = True ,margin = 50 ,
    #            bbox=(600,600),inline='None',autocurve = "True")
    
    
    
    # bbox = bbox.contract(20)
    
    myPlot.add(graph ,layout=layout,margin = 50,inline='None')
    myPlot.redraw()
    ctx = cairo.Context(myPlot.surface)
    
    ctx.set_font_size(25)
    
    # Choose the appropriate title
    if "label" in name_of_graph:
        title = "Sample: " + name_of_graph.split(sep="/")[0] + " - Label graph , " + "Gene: " + name_of_graph.split(" ")[0].split("_")[4]

    elif "threshold" in name_of_graph:
        
        title = "Sample: " + name_of_graph.split(sep="/")[0] + " - Threshold: " + name_of_graph.split(sep="_")[1] 
    elif "dynamic" in name_of_graph:
        title = "Sample: " + name_of_graph.split(sep="/")[0] + " - dynamic " 
    elif "normalized" in name_of_graph:
        title = "Sample: " + name_of_graph.split(sep="/")[0] + " - normalized"
    # elif "label" in name_of_graph:
    #     title = "Sample: " + name_of_graph.split(sep="/")[0] + " - Label graph , " + "Gene: " + name_of_graph.split(" ")[0].split("_")[4]

    drawer = TextDrawer(ctx, title, halign=TextDrawer.CENTER)
    drawer.draw_at(20, 25, width=550)
    #print("Saving image in : " + str(name_of_graph) + "\n")
    myPlot.save(name_of_graph)
    
    return 

def PlotDistribution(distance,sample):
    '''
        Input : similarity matrix
        -------------------------
        Output : figure object
    '''
    A = distance.flatten()
    A[A==0] = 1/10
    A = np.log10(A*10)
    # t = sns.distplot(A,hist=True, kde=True,
    #                  bins=int(300 / 5), color="darkblue",
    #                  hist_kws={'edgecolor': 'black'},
    #                  kde_kws={'linewidth': 2})
    t = sns.histplot(x=A, kde=True,stat="probability",
                     bins=int(300 / 5))
    #t.set_xlim(0,1)
    t.axes.set_title("Distribution of similarities - " + str(sample))
    figure = t.get_figure()
    return figure,A

def PlotGraphVisualStyle(clust,name_of_graph):   
    
    visual_style = dict()
    visual_style["bbox"] = (500,500)
    #visual_style["vertex_label"] = labels
    ig.plot(clust,name_of_graph,mark_groups = False,**visual_style)
    return


