'''
Created on May 21, 2014

@author: Timothy Ellersiek
'''

if __name__ == '__main__':
    pass

from dbLoader import loader
import plot as plot


graph = loader()
graph.edge_contraction()
graph.reduction_algorithm()


#plot.mpl_plot(graph)

#graph.edges_iter()
#graph.adj_iter()



  

    
