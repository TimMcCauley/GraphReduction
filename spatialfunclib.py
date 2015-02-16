'''
Created on Oct 7, 2014

@author: Timothy Ellersiek
'''

import math
import numpy as np
from scipy.spatial import distance
from itertools import combinations
from shapely.geometry import LineString
import collections
import ogr
import osr
from datetime import datetime
date_format = "%m/%d/%y %H:%M:%S"


#
# returns convex hull
#
def convex_hull(nodes):
       
    geomcol = ogr.Geometry(ogr.wkbGeometryCollection)
    
    for node in nodes:
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(float(node.coords[0]), float(node.coords[1]))
        
        source = osr.SpatialReference()
        source.ImportFromEPSG(4326)
        target = osr.SpatialReference()
        target.ImportFromEPSG(32632)
        transform = osr.CoordinateTransformation(source, target)
                        
        point.Transform(transform)
                
        geomcol.AddGeometry(point)
        
    convexHull = geomcol.ConvexHull()
    
    return convexHull.GetArea()
    
#
# returns the path bearing between two points specified in degrees. Formula adapted from: http://www.movable-type.co.uk/scripts/latlong.html
#
def path_bearing(a_lat, a_lon, b_lat, b_lon):
    
    a_lat = math.radians(a_lat)
    a_lon = math.radians(a_lon)
    b_lat = math.radians(b_lat)
    b_lon = math.radians(b_lon)
    
    y = math.sin(b_lon-a_lon) * math.cos(b_lat)
    x = math.cos(a_lat) * math.sin(b_lat) - math.sin(a_lat) * math.cos(b_lat) * math.cos(b_lon-a_lon)
    
    bearing = math.atan2(y, x)
    
    return math.fmod(math.degrees(bearing) + 360.0, 360.0)
      

#
# these two function checks if interesections between linestrings occur
#
def ccw(A,B,C):
    
    return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
    
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)


#
# returns bounding box of face
#
def bounding_box(face_nodes):
    ''' returns bounding box of polygon '''
    
    min_x = [(float(p.coords[0])) for p in face_nodes]
    min_x = min(min_x)
    
    max_x = [(float(p.coords[0])) for p in face_nodes]
    max_x = max(max_x)
    
    min_y = [(float(p.coords[1])) for p in face_nodes]
    min_y = min(min_y)
    
    max_y = [(float(p.coords[1])) for p in face_nodes]
    max_y = max(max_y)

    return [min_x,max_x,min_y,max_y]


#
# returns all permutations of given linestrings
#
def ls_permutations(all_linestrings):
    ''' return a list of all combinations
    of linestrings '''
    
    return list(combinations(all_linestrings, r=2))


#
# returns True if linestrings intersect
#
def intersect_all(perms):
    ''' intersect all projected linestrings with
    each other and return True if they do '''
    
    for i in range(len(perms)):
                    
        a_edge = LineString([perms[i][0][0],perms[i][0][1]])
        b_edge = LineString([perms[i][1][0],perms[i][1][1]])

        if a_edge.intersects(b_edge) == True:
            if a_edge.touches(b_edge) == False:
            
                print 'intersection'
              
                return True


#
# Calculates Hausdorff distance of polygon
# returns True if distance greater than hd_max
#
def hausDorffDistance(graph, s, e, face_nodes, hd_max):
    
    ''' creates two paths of nodes.. etc '''
    
    a = datetime.now()
    edgelist_a = collections.deque(face_nodes)
    edgelist_a.rotate(-face_nodes.index(s))
    edgelist_1 = list(edgelist_a)
    edgelist_1 = edgelist_1[edgelist_1.index(s):edgelist_1.index(e)+1]

    edgelist_b = collections.deque(face_nodes)
    edgelist_b.rotate(-face_nodes.index(e))
    edgelist_2 = list(edgelist_b)
    edgelist_2 = edgelist_2[edgelist_2.index(e):edgelist_2.index(s)+1]

    if len(edgelist_1) > len(edgelist_2):
        long = edgelist_1
        short = edgelist_2
    else:
        long = edgelist_2
        short = edgelist_1
    
    shortest_dist = []
    for node in long:
        
        source = osr.SpatialReference()
        source.ImportFromEPSG(4326)
        # target EPS
        target = osr.SpatialReference()
        target.ImportFromEPSG(32632)
        transform = osr.CoordinateTransformation(source, target)
        
        point = ogr.Geometry(ogr.wkbPoint)
        
        # this node coords
        this_coords = return_coords_from_graph(graph, node)
                
        point.AddPoint(float(this_coords[0]),float(this_coords[1]))
        point.Transform(transform)
        
        this_coords = [point.GetPoint()[0],point.GetPoint()[1]]
        
        coords_short = []
        
        for node_1 in short:
            point = ogr.Geometry(ogr.wkbPoint)
            coords = return_coords_from_graph(graph,node_1)
            point.AddPoint(float(coords[0]), float(coords[1]))
            point.Transform(transform)
            
            coords = [point.GetPoint()[0],point.GetPoint()[1]]
            
            coords_short.append(coords)
            
        distances_to_opposite = distance.cdist([this_coords],coords_short, 'euclidean')
                
        minval = np.min(distances_to_opposite)
    
        shortest_dist.append(minval)
        
        if minval > hd_max:
            b = datetime.now()
            delta = b - a
            
            return True
        
    b = datetime.now()
    delta = b - a
        
    #print 'HD distance meter ', np.max(shortest_dist)
    return False
    
# helper function to return coords from graph    
def return_coords_from_graph(graph,node):
        ''' helper function to get out_node coords 
        takes either node object or id as int '''
        
        try: 
            return graph.nodeID_index[node].coords
        
        except (AttributeError):
            
            return graph.nodeID_index[node.vID].coords
    

