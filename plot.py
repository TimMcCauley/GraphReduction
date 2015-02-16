'''
Created on May 21, 2014

@author: Timothy Ellersiek
'''

import matplotlib.pyplot as plt

def mpl_plot(graph, cnt=None, boundingbox=[], skel=[], longest_path=[]):
    
    #print 'skel ', skel  
    #print boundingbox
    
    #fig = plt.figure(figsize=(20, 8))
    fig = plt.figure()
    
    ax = plt.subplot(1,1,1)
    
    for k in graph.node_index.keys():    
        
        ax.scatter(k.coords[0], k.coords[1], s=8, marker='o', alpha=0.3)
        #ax.annotate(k.vID, xy=(k.coords[0], k.coords[1]), fontsize=12)

        for edge in graph.node_index[k]:
            #print 'edgeID ', edge.eID
            for node in graph.edge_index[edge]:
                #print 'nodeID ', node.vID
                if node.vID != k.vID:
                    
                    #print k.coords
                    #print node.coords
                    
                    #l, = ax.plot([k.coords[0],node.coords[0]], [k.coords[1],node.coords[1]], marker='o', alpha=0.1)
                    l, = ax.plot([k.coords[0],node.coords[0]], [k.coords[1],node.coords[1]], alpha=0.5)
                 
                    #string = [u.eID for u,g in graph.edge_index.iteritems() if u == edge][0]
                    #pos = [(k.coords[0]+node.coords[0])/2, (k.coords[1]+node.coords[1])/2]
                    #plt.text(pos[0], pos[1], string, size=10, rotation=12, color = l.get_color(), ha="center", va="center")
                     
                     
                    #plt.text(pos[0], pos[1], string, size=9, rotation=12, color = l.get_color(), #ha="center", va="center",bbox = dict(ec='1',fc='1'))
                    #ax.annotate(k.vID, xy=(k.coords[0], k.coords[1]), fontsize=12)
                    #ax.annotate(node.vID, xy=(node.coords[0], node.coords[1]), fontsize=12)

    
    if len(skel) > 0:  
        for line in skel:
            ax.plot([line[0][0],line[1][0]],[line[0][1],line[1][1]],'--', color='green', dashes=(1,2))
    
    if len(longest_path) > 0:
        for i in range(len(longest_path)-1):
            ax.plot([longest_path[i][1][0],longest_path[i+1][1][0]],[longest_path[i][1][1],longest_path[i+1][1][1]],color='orange')
        
    
    if len(boundingbox) > 0:
        ax.set_xlim(boundingbox[0],boundingbox[1])
        ax.set_ylim(boundingbox[2],boundingbox[3])
    
    
    if cnt != None:
        #fig.savefig('results_screens/'+'iteration'+str(cnt)+'.png',dpi=75)
       
        plt.show()
        #ax.axis('equal')
        plt.close(fig)
    
    else: 
        #ax.axis('equal')
        plt.show()
    
    
    
    