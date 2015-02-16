'''
Created on May 21, 2014

@author: Timothy Ellersiek
'''

import os
import itertools
import ogr
import osr

#===============================================================================
# Node
#===============================================================================
class Node:
    """
    Node of the graph
    3D point + list of neighboring arcs
    """
    
    count = 0
    
    def __init__(self,p):
        self.value      = p
        #print 'self value ', self.value

        self.id         = Node.count
        Node.count      = Node.count + 1

        self.neighbors  = []
        self.visited    = False
        
    def __eq__(self,other):
        return self.value == other.value

    def __repr__(self):
        return self.id.__repr__()

            
    def add_neighbor(self,a):
        try :
            return self.neighbors.index(a)
        except:
            self.neighbors.append(a)
            return len(self.neighbors)-1
    

#===============================================================================
# Arc
#===============================================================================
class Arc:
    """
    An arc between two nodes
    """
    
    def __init__(self,p1,p2):
        self.p1 = p1
        self.p2 = p2
        
    
    def __eq__(self,other):
        """
        Compare two arcs
        As the arcs represent undirected links, arc(p1,p2) = arc(p2,p1)
        """
        return ((self.p1 == other.p1) and (self.p2 == other.p2)) or ((self.p1 == other.p2) and (self.p2 == other.p1))
        
    def __repr__(self):
        return self.p1.__repr__() +" -- "+ self.p2.__repr__()
    
    def get_other_end(self,p):
        if p==self.p1:
            return self.p2
        elif p==self.p2:
            return self.p1
        else:
            print 'O_o'


#===============================================================================
# Graph_proj
#===============================================================================
class Graph:
    """
    This class provide tools to build a undirected graph, to manipulate it and
    to find minimal path
    """
        
    def __init__(self,nodes=[],arcs=[]):
        Node.count  = 0
        self.nodes  = nodes
        self.arcs   = arcs
        self.adj = {}
        
    def reset_visit(self):
        for node in self.nodes:
            node.visited = False
    
    def find(self,id):
        """ Find the node of the graph which value is 'value' """
        result = False
        for node in self.nodes:
            if node.id == id:
                result = node
                break
        return result
        
    #===========================================================================
    # Change the topology of the graph
    #===========================================================================
    def add_arc(self,p1,p2):
        """
        @arg p1 : Node
        @arg p2 : Node
        @return : the created arc itself
        """
        
        a = Arc(p1,p2)
        if p1==p2:
            return None
        
        if not a in self.arcs:
            self.arcs.append(a)    
            p1.add_neighbor(a)
            p2.add_neighbor(a)
        
        # filling up adjacency list
        if not p1.id in self.adj:
            self.adj[p1.id] = [p2.id]
        else:
            self.adj[p1.id].append(p2.id)
        
        if not p2.id in self.adj:
            self.adj[p2.id] = [p1.id]
        else:
            self.adj[p2.id].append(p1.id)
        
        return a
        
    def add_node(self,p):
        """
        @arg p : the value of the node to be added. The information contained in the node
                 to be created. (point, index, whatever)
        
        @return : the node that had been created, or the one that already
                  existed with the given value p
        """
        node_p = Node(p)
        try :
            index = self.nodes.index(node_p)
            Node.count -= 1 #take care, by creating a node, we have increment the count of existing Nodes.
            return self.nodes[index]
        except:
            self.nodes.append(node_p)
            return node_p
        
        
    def clean_skeleton(self):
        
        self.nodes = []
        self.arcs = []
        self.adj = {}
        Node.count  = 0
            
        
    #mine
    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        """ find all paths from start_vertex to 
            end_vertex in graph """
        graph = self.adj
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            #print [path] 
            return [path]
        if start_vertex not in graph:
            return []
        paths = []
        
        for vertex in graph[start_vertex]:
            
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex, 
                                                     end_vertex, 
                                                     path)
                
                for p in extended_paths: 
                    paths.append(p)  
        return paths
    
    
    def get_skeleton_outer_nodes(self):
        """ return all points of degree one in 
        skeleton, these are the outest points """
        
        outer_nodes = []
        for k,v in self.adj.iteritems():
            if len(v) == 1:
                outer_nodes.append(k)
        
        return outer_nodes
    
    
    def find_longest_path(self):
        ''' finds longest path from degree 1 to degree 1 nodes
        in skeleton and returns IDs and coordinates '''
        
        deg1_nodes = self.get_skeleton_outer_nodes()
        
        permutations = list((set(s) for s in itertools.permutations(deg1_nodes, 2) if len(set(map(abs, s))) == len(s)))
            
        new_perms = []
        for perm in permutations:
            if perm not in new_perms:
                new_perms.append(perm)
                
        permutations = [list(perm) for perm in new_perms]
        
        # transform WGS to UTM
        source = osr.SpatialReference()
        source.ImportFromEPSG(4326)
        target = osr.SpatialReference()
        target.ImportFromEPSG(32632)
        transform = osr.CoordinateTransformation(source, target)
        
        max_distance = 0
        
        # all permutations
        for tuple in permutations:  
            paths = self.find_all_paths(tuple[0], tuple[1])
            line = ogr.Geometry(ogr.wkbLineString)
            for path in paths:
                for i in range(len(path)):
                    this_node = self.find(path[i])
                    line.AddPoint(float(this_node.value[0]),float(this_node.value[1]))
                line.Transform(transform)
                #print "Length = %d" % line.Length()
                if line.Length() > max_distance:
                    max_distance = line.Length()
                    longest_path = path
        longest_path_coords = [] 
        for node in longest_path:
            longest_path_coords.append((self.find(node).id,self.find(node).value))
            
        return longest_path_coords
    
        
    def find_1deg(self,start_node, path = [], medial_node = None):
        ''' recursive function to find
        all 1 degree nodes '''

        paths = []
    
        # add starting node which lies on longest_path of skeleton
        if medial_node is not None:
            path.append(medial_node)
    
        # if this node has degree 3 then add to path
        if len(self.adj[start_node]) == 3:
            path = path + [start_node]
    
            # loop through nbr nodes of this node
            for nbr in self.adj[start_node]:
    
                # if not in path
                if nbr not in path:
    
                    # start algorithm again
                    extended_paths = self.find_1deg(nbr, path)
                    for p in extended_paths: 
                        
                        paths.append((p[0], p[-1]))
        # else dead end is found = outer node on skeleton graph              
        else:
            path = path + [start_node]
            paths.append((path[0], path[-1]))
            
        return paths
    
    
    def projection_tuples(self,longest_path):
        ''' in combination with find_1deg() this function
        returns all 1 degree nodes connected to
        longest path of skeleton '''
        
        projection_tuples = []
        #grab ids only
        longest_path = [step[0] for step in longest_path]
                
        for step in longest_path[1:-1]:
            #print 'node ', step
            if len(self.adj[step]) == 3:
                for nbr in self.adj[step]:
                    if not nbr in longest_path:
                        projection_tuples.append(self.find_1deg(nbr,[],step))

        return projection_tuples


    def get_coords(self,tuples_list,r):
        ''' grab coordinates of degree 1 skeleton nodes 
        and add them to index '''
        
        medial_1dg_idx = {}
        for r_id, step in zip(r,tuples_list):
            for tuple in step:
                
                if (r_id,tuple[0]) in medial_1dg_idx:
                    medial_1dg_idx[(r_id,tuple[0])].append(self.find(tuple[1]).value)
                else:    
                    medial_1dg_idx[(r_id,tuple[0])] = [self.find(tuple[1]).value]
        
        return medial_1dg_idx
    
    
#===============================================================================
if __name__ == '__main__':
    g = Graph()

    linestrings = [[(1.696, 48.145), (1.81395837657777, 48.3794541577622)], [(2.067, 48.242), (1.81395837657777, 48.3794541577622)], [(1.911, 48.713), 
    (1.77037164893532, 48.4983758587174)], [(1.483, 48.674), (1.77037164893532, 48.4983758587174)], 
    [(1.81395837657777, 48.3794541577622), (1.77037164893532, 48.4983758587174)]]
    #test [(1.911, 48.713), (1.81395837657777, 48.3794541577622)]

    for i in range(len(linestrings)):
        p1 = g.add_node(linestrings[i][0])
        p2 = g.add_node(linestrings[i][1])
        g.add_arc(p1,p2)

    print 'arcs ', g.arcs
    print 'node count ', g.nodes
    print 'adja ', g.adj
    #for node in g.nodes:
        #print  node, node.value
        #print 'neighbors ', node.neighbors

    #print g.find_all_paths(5, 1)
      
    # get all skeleton outer points with length 1
    #deg1_nodes = g.get_skeleton_outer_nodes()
    
    #longest_path = g.find_longest_path()
    
    # grab outer nodes

    #paths_to_outer_nodes = g.paths_to_1g(longest_path[1:-1])
    
