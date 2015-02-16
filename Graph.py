'''
Created on May 21, 2014

@author: Timothy Ellersiek '''

import spatialfunclib as spatial
import bisect
import operator
import plot as plot
import numpy as np
import copy
import sys
import psycopg2 as pg
import ogr
import re
from decimal import *
from skeleton_graph import skel
from collections import defaultdict
import os
import cProfile

class node():
    ''' Node object of graph '''
    
    def __init__(self, vID, coordinates):
        self.vID = vID
        self.coords = coordinates

    def __hash__(self):
        return hash(self.vID)

    def __eq__(self, other):
        return self.vID == other.vID
    
     
class edge():
    ''' Edge object of graph '''
    
    def __init__(self, nodeList, eID):
        self.nodes = nodeList
        self.eID = eID
    
    def __hash__(self):
        return hash(self.eID)

    def __eq__(self, other):
        return self.eID == other.eID
    
    
class multiGraph():
    ''' Graph object containing all steps for reduction '''
    
    def __init__(self):
        self.e_cnt = 0
        self.n_maxID = None
        
        self.node_index = {}   # empty adj list
        self.edge_index = {}    # empty edge index
        self.nodeID_index = {} # emptry nodeid index

    def add_edge(self,x,y,x_coords=[],y_coords=[]):
        ''' Adding an edge to the graph:
        a) creates an edge object and adds the edge id (incremental) 
        to edge_index (keys) with its two containing nodes (values)
        b) creates two node objects for in- and out-node which 
        are added to node_index dict -> key is the node object and 
        value is the edge it belongs to. This executes twice as edge
        consists of two nodes. New edges are inserted into the value list
        according to their bearing in and increasing order (clockwise), e.g.
        
        node_index dictionary structure:
        NODE (objects) | EDGES (objects)
        1              | 1
        2              | 1
        
        edge_index dictionary structure:
        EDGE (edge id) | NODES (node ids as list)
        1              | 1,2 '''
        
        # x and y can not be the same
        #print x,y
        if x == y:
            sys.exit()
                             
        # create an edge object
        this_edge = edge(list([x,y]),self.e_cnt)
        
        # create an x node object (in node)
        x_node = node(x,x_coords)
        
        # check if x_node is in node_index, if it is insert to edge to list maintaining clockwise order
        if x_node in self.node_index: 
            # calculate bearing of this edge
            this_projection_bearing = spatial.path_bearing(x_coords[1], x_coords[0], y_coords[1], y_coords[0])
                        
            # calculate all bearings of other edges connected to this node
            other_bearings_keys = self.calc_bearings(x_node)
            
            # this is the edge which has to be inserted
            insert_new_edge_with_bearing = [this_edge,this_projection_bearing]
            
            # use bisect library to insert according to bearing          
            self.node_index[x_node].insert(bisect.bisect_left(other_bearings_keys,insert_new_edge_with_bearing[1]),insert_new_edge_with_bearing[0])

        # if x_node is not in node_index, add it
        else:
            self.node_index[x_node] = [this_edge]
            self.nodeID_index[x] = x_node
            
        # and now repeat for y node - create an y node object (out node)
        y_node = node(y,y_coords)
        # check if y_node is in node_index, if it is insert to edge to list maintaining clockwise order
        if y_node in self.node_index:
            
            # calculate bearing of this edge
            this_projection_bearing = spatial.path_bearing(y_coords[1], y_coords[0], x_coords[1], x_coords[0])
            
            # calculate all bearings of other edges connected to this node
            other_bearings_keys = self.calc_bearings(y_node)

            # this is the edge which has to be inserted
            insert_new_edge_with_bearing = [this_edge,this_projection_bearing]

            # use bisect library to insert according to bearing          
            self.node_index[y_node].insert(bisect.bisect_left(other_bearings_keys,insert_new_edge_with_bearing[1]),insert_new_edge_with_bearing[0])

        # if y_node is not in node_index, add it
        else:
            self.node_index[y_node] = [this_edge]
            self.nodeID_index[y] = y_node
    
        # add this edge object to edge_index
        self.edge_index[this_edge] = [x_node, y_node]
        
        # increment edge_counter        
        self.e_cnt += 1
        
        return
    
    
    def edge_contraction(self):
        ''' removes all vertices of degree 1 recursively '''
        
        again = False
        
        for k in self.node_index.keys():
            if len(self.node_index[k]) == 1:
                again = True
                # remove edge from other node of edge
                for node in self.node_index[k][0].nodes:
                    if node != k.vID:
                        for t,s in self.node_index.iteritems():
                            if t.vID == node:
                                self.node_index[t].remove(self.node_index[k][0])
                                
                # edge removal from edge_index           
                del self.edge_index[self.node_index[k][0]]
                del self.node_index[k]            
    
        if again == True: 
            self.edge_contraction()
        else:
            print 'Edges contracted..finished'
            return
            
    
    def reduction_algorithm(self):
        ''' builds connection to Straight Skeleton database 
        starts algoritm until no faces are found (condition = false) '''
        
        try:
            conn = pg.connect("dbname='DB_0511' user='sftest' host='localhost' password='sftest'")
        except:
            print "Unable to connect to PostGIS CGAL the database"
          
            
        condition = True
        started_cnt = 0
        
        while condition:
            
            print '\nStarting Morris V2 for the', started_cnt, 'time'
            condition = self.find_faces(conn, started_cnt)
                            
            if condition == False:
                
                print '\nMorris V2 completed'
                
                # save poly lines at the end
                nodes_cp = copy.deepcopy(self.node_index)
                edges_cp = copy.deepcopy(self.edge_index)
                self.save_polylines(0,nodes_cp, edges_cp)
                
                conn.close()
                break
            
            started_cnt += 1
                          
                       
    def find_faces(self, conn, started_cnt):  
        ''' iterates all vertices and its edges
        and finds closed faces '''
        
        cycle = False
        cnt = 0 
        reduced_in_this_iteration = False
        
        for k in self.node_index.keys():
            if cnt % 1000 == 0:
                print cnt,'/',len(self.node_index)
            
                self.save_to_csv(started_cnt,cnt)
            
            if k in self.node_index.keys():
                
                # nodes with degree 3 or larger as starting points
                # to find faces
                if len(self.node_index[k]) >= 3:
                    
                    # function to iterate edges of this node
                    def iterateEdges(cycle,cnt,started_cnt):
                        
                        faceReducedBool = False
                            
                        # iterate through ordered edge list connected to vertice
                        for edge in self.node_index[k]:
                            
                            face_edge_objects = []
                            face_node_objects = []
                            
                            next_edge = edge
                            in_node = k
                            
                            if len(self.node_index[in_node]) == 1:
                                print 'Dead End, breaking..'
                                break
                            
                            face_edge_objects.append(edge)
                            face_node_objects.append(in_node)
                            
                            # finding face cycle
                            while not cycle:
                                
                                for next_node in self.edge_index[next_edge]:
                                    if next_node.vID != in_node.vID:
                                        out_node = next_node
                                                                                            
                                next_edge_list = self.node_index[out_node]
                                                                
                                # dead end?
                                if len(next_edge_list) == 1:
                                    print 'Dead End, breaking..'
                                    break
                                                                
                                # get index of this edge in edges_list of node
                                for e in next_edge_list:
                                    if e.eID == next_edge.eID:
                                        index_origin = next_edge_list.index(e)
                                                   
                                next_edge = self.rightmost_edge(next_edge_list, index_origin)
                                
                                # zipfel may occur
                                check1 = self.nodeID_index[next_edge.nodes[0]]
                                check2 = self.nodeID_index[next_edge.nodes[1]]
                                if len(self.node_index[check1]) == 1:
                                    break
                                if len(self.node_index[check2]) == 1:
                                    break
                                
                                in_node = out_node
                                
                                # check if face is found
                                cycle = self.evaluate_end_condition(k.vID, out_node.vID)
                                    
                                if cycle is True:
                                    
#                                     print 'cycle found ', [nodeness.vID for nodeness in face_node_objects]
                                                                
                                    # if a node occurs twice in the face then return False 
                                    if len(set([nodeness.vID for nodeness in face_node_objects])) < len([nodeness.vID for nodeness in face_node_objects]):
                                        print 'Duplicates in node list of face, breaking..'
                                        cycle = False
                                        break
                                    
                                    # try to reduce face                   
                                    reduced = self._reduction(face_edge_objects, face_node_objects, conn)                                             
                                    
                                    if reduced == False:
                                        cycle = False
                                        break
                                        #return
                                        
                                    # if reduced is True, set cycle false and 
                                    # return from iterateEdges and continue to next node
                                    else: 
                                        cycle = False
                                        faceReducedBool = True
                                        break
                                        #return faceReducedBool
                                    
                                face_edge_objects.append(next_edge)
                                face_node_objects.append(in_node)
                    
                        # if no face is found within the edges of this node return  
                        return faceReducedBool
                     
                    # cProfile
#                     def wrapper(b):
#                         b.append(iterateEdges(cycle,cnt,reduced_in_this_iteration))
#                      
#                     b = []
#                     cProfile.runctx('wrapper(b)', globals(), locals(), sort=1)
#                     reduced_in_this_iteration = b[0]

                    face_has_been_reduced = iterateEdges(cycle,cnt,started_cnt)
                    
                    if face_has_been_reduced == True:
                        reduced_in_this_iteration = True
                                                 
            cnt += 1
                    
        # if any face was reduced in this iteration then return True and restart algo
        if reduced_in_this_iteration == True:
            print 'node_index end reached, restarting algorithm..'
            return True
        
        # if no cycles are found after iterating the entire node_index, 
        # then return False, algorithm is completed
        if cycle == False:
                                                            
            self.save_to_csv(started_cnt,None)
    
            return False   
        
           
    def evaluate_end_condition(self,vID, out_node):
        ''' if out_node is start_node then cycle is found '''
        
        if out_node == vID:
            return True
        else:
            return False
    
    
    def save_to_csv(self,started_cnt,cnt=None):
        ''' saves edge list to csv as WKT '''
        
        if cnt == None:
            ncnt = str('_fin')
            
            #plot.mpl_plot(self,cnt)
            
            self.save(started_cnt,ncnt)
            
        else:
            ncnt = str(cnt)
        
            if cnt % 1000 == 0:
                
                #plot.mpl_plot(self,cnt)
                
                self.save(started_cnt,ncnt)
                                     
    
    def save(self,started_cnt,ncnt=None):
        
        filename = str(started_cnt)+'_'+'results_morrisV2_'+ncnt+'.csv'
        f = open(os.path.join('/home/giscience/workspace/Morris2014/results/paper', filename), 'w')
        
        for k,v in self.edge_index.iteritems():
            linestring = str('LINESTRING(')
            for i in range(len(v)):
                for latlon in v[i].coords:
                    linestring += str(latlon)
                    linestring += str(' ')
                if i == 0:
                    linestring += str(',')
            linestring += str(')')
            linestring += str(';')
            linestring += str(v[0].vID)
            linestring += str(';')
            linestring += str(v[1].vID)
            linestring += str('\n')
            f.write(linestring)
                                  
        f.close()
        
          
    def save_polylines(self, counter=None, nodes_cp=[], edges_cp=[]):
        ''' function cuts the graph up into linestrings 
        between intersections finds degree 3 or larger 
        nodes as starting points for linestrings, if none 
        are found then degree 1 as starting point '''
        
        nodes_of_linestring = []
        
        in_node = None 
        deg3 = False
        
        for k,v in nodes_cp.iteritems():
            if len(nodes_cp[k]) >= 3:                                                                                      
                in_node = k
                deg3 = True
                break
        
        # get last line
        if deg3 == False:
            for c,f in nodes_cp.iteritems():
                if len(nodes_cp[c]) == 1:                                                                                      
                    in_node = c
                    break
        
        if in_node != None:
                                   
            next_edge = nodes_cp[in_node][0]
            nodes_of_linestring.append(in_node.coords)  
                           
            while 1:
                
                del_nodes = []
                nodes_cp[in_node].remove(next_edge)
                                                                  
                for next_node in edges_cp[next_edge]:
                    if next_node.vID != in_node.vID:
                        out_node = next_node
                
                # Polyline completed? Intersection found?
                if len(nodes_cp[out_node]) >= 3:
                                        
                    nodes_of_linestring.append(out_node.coords)
                    nodes_cp[out_node].remove(next_edge)
                    
                    for dn in del_nodes:
                        nodes_cp.pop(dn, None)
                        
                    filename = 'polyline_results_'+str(counter)+'.csv'
                    f = open(os.path.join('/home/giscience/workspace/Morris2014/results/LineStrings', filename), 'w')
                    linestring = str('LINESTRING(') 
                    for coord in nodes_of_linestring:
                        linestring += str(coord[0])
                        linestring += str(' ')
                        linestring += str(coord[1])
                        linestring += str(',')
                    linestring = linestring[0:-1]
                    linestring += str(')')                           
                    f.write(linestring)
                    f.close()
                    
                    counter += 1
                    
                    edges_cp.pop(next_edge, None)
                    break  
                
                # Or just a dead end?
                if len(nodes_cp[out_node]) == 1:
                    
                    nodes_of_linestring.append(out_node.coords)
                    nodes_cp[out_node].remove(next_edge)
                    
                    del_nodes.append(out_node)
                    for dn in del_nodes:
                        nodes_cp.pop(dn, None)
                    
                    filename = 'polyline_results_'+str(counter)+'.csv'
                    f = open(os.path.join('/home/giscience/workspace/Morris2014/results/LineStrings', filename), 'w')
                    linestring = str('LINESTRING(')
                    for coord in nodes_of_linestring:
                        linestring += str(coord[0])
                        linestring += str(' ')
                        linestring += str(coord[1])
                        linestring += str(',')
                    linestring = linestring[0:-1]
                    linestring += str(')')                           
                    f.write(linestring)
                    f.close()    
                                
                    counter += 1
                    
                    edges_cp.pop(next_edge, None)
                    break 
                
                # remove edge from edges                   
                edges_cp.pop(next_edge, None)
                # append out_node coords to linestring
                nodes_of_linestring.append(out_node.coords)
                # add out_node to list of nodes to be deleted
                del_nodes.append(out_node)
                # remove next_edge from out_node node
                nodes_cp[out_node].remove(next_edge)
                in_node = out_node
                # get next edge
                next_edge = nodes_cp[out_node][0]
            
            # start function again
            self.save_polylines(counter, nodes_cp, edges_cp)             
                
    
    def _remove_edges(self,edges, nodes):
        # remove edges of face after reduction, each edge occurs twice as they are
        # connected to two nodes, remove from this node and last node from node_index
        # also remove edge from edge_index
        for i,node in enumerate(nodes):
            
            idx_edge = edges.index(edges[i])
            
            #print 'removing edges in remove_edges().. ', edges[i].eID, edges[(idx_edge + len(edges)-1)%len(edges)].eID
            # 1st removal node_index
            self.node_index[node].remove(edges[i])
            # 2nd removal node_index
            self.node_index[node].remove(edges[(idx_edge + len(edges)-1)%len(edges)])
            
            # edge removal from edge_index           
            del self.edge_index[edges[i]]
            
        
    def rightmost_edge(self,ordered_edges_list, idx):
        ''' returns rightmost edge '''
        
        return ordered_edges_list[(idx + len(ordered_edges_list)-1)%len(ordered_edges_list)] 
        
    
    def calc_bearings(self,this_node):
        ''' calculates bearings of all 
        other edges connected to this_node '''
        
        brngs = []
        # iterate through all edge objects of this_node in node_index
        for check_edge in self.node_index[this_node]:
            
            # for the bearing we need to know the direction, as we are always starting from in_node or out_node
            if this_node.vID != check_edge.nodes[0]:
                check_edge.nodes.reverse()
             
            # create in_node_object and out_node_object
            #in_node_object_old = [k for k,v in self.node_index.iteritems() if k.vID == check_edge.nodes[0]][0]
            in_node_object = self.nodeID_index[check_edge.nodes[0]]
            
            #out_node_object_old = [k for k,v in self.node_index.iteritems() if k.vID == check_edge.nodes[1]][0]
            out_node_object = self.nodeID_index[check_edge.nodes[1]]
            
            # calculate bearing
            compare_projection_bearing = spatial.path_bearing(in_node_object.coords[1],in_node_object.coords[0],out_node_object.coords[1], out_node_object.coords[0])
            
            # add bearing to bearing list
            brngs.append(compare_projection_bearing)
           
        return brngs
    
        
    def increment_max_nodeID(self):
        ''' initiate new node ID '''

        return (max(k.vID for k,v in self.node_index.iteritems()))+1
        
        
    def _reduction(self,face_edges_objs, face_nodes_objs, conn):
        ''' this reduction function contains several steps 
        at first it checks if the convex hull exceeds a given
        threshold, returns False if it does.
        Then it fires a sql statement to the CGAL library to return
        the straight skeleton of the polygon. Afterwards 
        start and end point of the polygon are found 
        which represent the longest path through it,
        with these one may calculate the Hausdorff distance
        which however is not needed if you use the convex hull.
        Then the reduced polyline R is added to the polygon.
        From R (which is part of the skeleton) we can derive the endpoints
        of the skeleton graph (degree 1) and compare these to the polygon 
        coordinates to know which have to be connected to R. '''
        
        # build convex hull, return area and if bigger than given theshold return false 
        # (10000 for road network graph, 2000 without road network)
        convexHullThreshold = 10000
        if spatial.convex_hull(face_nodes_objs) > convexHullThreshold:
            print 'hull greater ', convexHullThreshold, ' square meters, return False'
            return False
         
        # build polygonWkt as String concatenation for ST_StraightSkeleton function
        polygon = self.buildWktPolygon(face_nodes_objs)                    
        
        # calculates bounding box of polygon     
#         bb = spatial.bounding_box(face_nodes_objs)
        #plots graph and zooms to bounding box
#         plot.mpl_plot(self,bb)
                                  
        # query db for skel here, return False if no skeleton can be calculated, 
        # meaning that face will be simply ignored
        cur = conn.cursor()
        try:
            cur.execute("""SELECT ST_AsText((ST_Dump(ST_StraightSkeleton(ST_PolygonFromText(%s)))).geom);""", (polygon,))
        except pg.Error as e:
            print e
            #print polygon
            #print [edgeness.eID for edgeness in face_edges_objs]
            #print [nodeness.vID for nodeness in face_nodes_objs]
            conn.rollback()
            #return True
            return False
        
        # fetch line coordinates of skel to build skel graph
        skeleton_linestrings = self.grabLineCoords(cur)        
        
        # returns all permutations of linestrings in skeleton
#         perms_skeleton_ls = spatial.ls_permutations(skeleton_linestrings)        
#         if spatial.intersect_all(perms_skeleton_ls) == True:
#             return False
        
        # create skeleton Graph_proj here and empty it if full
        sg = skel.Graph()
        sg.clean_skeleton()
        
        for i in range(len(skeleton_linestrings)):
            p1 = sg.add_node(skeleton_linestrings[i][0])
            p2 = sg.add_node(skeleton_linestrings[i][1])
            sg.add_arc(p1,p2)

        # find the longest path in skeleton graph
        longest_path = sg.find_longest_path()
        
        # match first and last coords of longest_path to 
        # coords of polygon to find start and end node    
        start, end = self.startEndIDsLongestPath(longest_path, face_nodes_objs)
        
#         if spatial.hausDorffDistance(self, start, end, [node.vID for node in face_nodes_objs],300) == True:
#             print 'HausDorff distance greater 300 meters, returning'
#             return False
        
        # plot graph with bounding box, skeleton and longest path of skeleton
#         plot.mpl_plot(self,None,spatial.bounding_box(face_nodes_objs),skeleton_linestrings,longest_path)

        # find projection tuples and match decimals to polygon decimals
        proj_tuples = sg.projection_tuples(longest_path)
                
        #add R
        R = self.add_R(start, end, longest_path[1:-1])
        
        # grab coordinates of degree 1 nodes connected to medial axis of skeleton
        proj_coordinates = sg.get_coords(proj_tuples,R[1:-1])     
                        
        # match coordinates of 1g nodes of skeleton to polygon
        matched_coordinates = self.match_1dg_skel_to_poly(face_nodes_objs, proj_coordinates)

        # remove edges from polygon
        self._remove_edges(face_edges_objs, face_nodes_objs)
        
        # check if to project nodes have nbrs, if yes : project, if not: delete node
        self.project_by_connect(matched_coordinates)
        
        # clean up dead ends of face
        self.clean_dead_ends(R)
          
        # face was reduced, return True  
        return True 
      

    def clean_dead_ends(self,R,cnt=None):  
        ''' cleans dead ends from skeleton which has been added to graph
        calls itself again to start from the other end of skeleton linestring '''
        
        for step in R:
            nodeObject = self.nodeID_index[step]
            if nodeObject in self.node_index:
                if len(self.node_index[nodeObject]) == 1:
                    # remove edge from other node of edge
                    for node in self.node_index[nodeObject][0].nodes:
                        if node != step:
                            key = self.nodeID_index[node]
                            self.node_index[key].remove(self.node_index[nodeObject][0])
                            #for t,s in self.node_index.iteritems():
                                #if t.vID == node:
                                    #self.node_index[t].remove(self.node_index[k][0])
                    # edge removal from edge_index           
                    del self.edge_index[self.node_index[nodeObject][0]]
                    del self.node_index[nodeObject]
                elif cnt == 1:
                    return
                else: 
                    R.reverse()
                    self.clean_dead_ends(R,1)
                                    
            
    def get_nbr_edges(self,this_node):
        ''' helper function to get nbr edges 
        takes either node object or id as int '''       
        
        try: 
            return self.node_index[this_node]
        except (AttributeError):
            key = self.nodeID_index[node]
            return self.node_index[key] 

    
    def get_out_node_coords(self,node_to_project):
        ''' helper function to get out_node coords 
        takes either node object or id as int '''
        
        return self.nodeID_index[node_to_project].coords
                            
        
    def get_out_node(self,in_node_id,this_edge):
        ''' helper function to get out_node object 
        takes either node object or id as int '''
        
        for nbr_node in self.edge_index.get(this_edge):
            try:
                if nbr_node.vID != in_node_id.vID:
                    return nbr_node
            except (AttributeError):
                if nbr_node.vID != in_node_id:
                    return nbr_node
                        
    
    def project_by_connect(self, list_to_project):
        ''' 1. adds edges from nodes with degree greater 
        2 to medial axis
        2. removes old edges afterwards '''
                
        for k,v in list_to_project.iteritems():
            medial_coords = self.get_out_node_coords(k[0])
            for proj_node in v:
                # get neighbors
                nbr_edges = copy.deepcopy(self.node_index[proj_node])
                # check if to project node has neighbors
                if len(nbr_edges) > 0:
                    self.add_edge(k[0], proj_node.vID, medial_coords, proj_node.coords)
                else:
                    del self.node_index[proj_node]
                        
    
    def match_1dg_skel_to_poly(self, poly_nodes, proj_coords):
        ''' match degree 1 nodes coords with polygon
        to get corresponding node ids from polygon  '''
        
        degree1_matched = {}
        
        for k,v in proj_coords.iteritems():
            for coords in v:
                for poly_node in poly_nodes:
                    if poly_node.coords == coords:
                        if k in degree1_matched:
                            degree1_matched[k].append(poly_node)
                        else:
                            degree1_matched[k] = [poly_node]
                                                      
        return degree1_matched
                        
    
    def startEndIDsLongestPath(self,longest_path, face_nodes_objs):
        ''' Match Start and End Coordinates of Longest Line
        in Skeleton Graph_proj to Coordinates in Polygon to find
        matching points in Polygon '''
        
        startNode = None
        endNode = None
    
        for coordinate in [longest_path[0],longest_path[-1]]:
            for node in face_nodes_objs:
                if coordinate[1] == node.coords:
                    if startNode is None:
                        startNode = node.vID
                    else:
                        endNode = node.vID
        
        return startNode, endNode
    
    
    def grabLineCoords(self, cur):
        ''' Grab Line Coordinates of StraightPolygon from
        DB to Decimal '''
        
        linestrings = []
        #loop through cursor
        for record in cur:
            record_string=str(record)
            #get edge pairs (numbers start always at pos 11)
            split1 = record_string[11+2:record_string.find(")")]            
            #split edges to single vertices
            split2 = re.split(",", split1)
            #split split2[0] again to declare decimals
            split3 = split2[0].split(" ")
            a1 = Decimal(str(split3[0]))
            a2 = Decimal(str(split3[1]))
            #split split[1] again to declare decimals
            split4 = split2[1].split(" ")
            b1 = Decimal(str(split4[0]))
            b2 = Decimal(str(split4[1]))
             
            #insert into linestrings array
            linestrings.append([(a1,a2),(b1,b2)])
        
        return linestrings
    
    
    def buildWktPolygon(self, face_nodes_objs):     
        ''' Builds a Wkt String, doesn't work with ogr
        due to float point differences when passing
        polygon coordinates to ring '''
        
        polygon = str('POLYGON((')

        for i in range(len(face_nodes_objs)):
            #print face_nodes_objs[i].coords[0], face_nodes_objs[i].coords[1]
            polygon+= str(face_nodes_objs[i].coords[0])
            polygon+= str(' ')
            polygon+= str(face_nodes_objs[i].coords[1])
            polygon+= str(',')
            if i==len(face_nodes_objs)-1:
                #print face_nodes_objs[0].coords[0], face_nodes_objs[0].coords[1]
                polygon+=str(face_nodes_objs[0].coords[0])
                polygon+=str(' ')
                polygon+=str(face_nodes_objs[0].coords[1])
                polygon+=str('))')
        
        return polygon
     
    
    def add_R(self, start, end, centroids):
        ''' Adds reduced line R from start to end point of face
        returns R '''  
        
        self.n_maxID += 1      
        R = []
        startCoords = self.nodeID_index[start].coords
        endCoords = self.nodeID_index[end].coords
        start = (start,startCoords)
        end = (end,endCoords)
        centroids = [start] + centroids + [end]
                
        for i in range(len(centroids)):
            # start 
            if i == 0:
                self.add_edge(centroids[i][0], self.n_maxID, centroids[i][1], centroids[i+1][1])
                R.append(centroids[i][0])
            # end
            elif i == len(centroids)-2:
                self.add_edge(self.n_maxID, centroids[i+1][0], centroids[i][1], centroids[i+1][1])
                R.append(self.n_maxID)
                R.append(centroids[i+1][0])
            # inbetween
            elif i < len(centroids)-2:
                self.add_edge(self.n_maxID, self.n_maxID+1, centroids[i][1], centroids[i+1][1])
                R.append(self.n_maxID )
                self.n_maxID += 1
        
        return R
            
            
    def adj_iter(self):
        ''' iterates through node_index 
        and yields adjacency matrix '''
        
        for k,v in self.node_index.iteritems():
            print k.vID, k.coords, ': ', [edges.eID for edges in v]          
    
    
    def edges_iter(self):
        ''' iterates through edge_index 
        and yields in_node / out_node '''
        
        for k,v in self.edge_index.iteritems():
            print k.eID,v