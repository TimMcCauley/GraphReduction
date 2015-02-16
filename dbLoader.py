'''
Created on May 21, 2014

@author: Timothy Ellersiek

usage example:
mg = multiGraph() 
mg.add_edge(1, 8, [54.3, 2.8], [53.3, 2.3])
mg.add_edge(3, 4, [52.3, 2.6], [53.1, 2.2])
mg.add_edge(3, 5, [52.3, 2.6], [53.1, 2.2])
mg.add_edge(3, 6, [52.3, 2.6], [53.1, 2.2])
mg.add_edge(3, 7, [52.3, 2.6], [53.1, 2.2])
mg.add_edge(3, 8, [52.3, 2.6], [53.1, 2.2])
mg.add_edge(3, 8, [52.3, 2.8], [53.1, 2.2]) '''


import psycopg2 as pg
import Graph as g
import copy
from decimal import *


def loader():

    try:
        #conn = pg.connect("dbname='taplus_testdb_20140327' user='tellersiek' host='teleagroplus.geog.uni-heidelberg.de' password='kiesmit2014'")       
        #FROM  new_dataset_jo_single_ls_ss_ss_topo.edge_data a
        conn = pg.connect("dbname='TeleAgroSubsets' user='sftest' host='localhost' password='sftest'")       

    except:
        print "I am unable to connect to the database"
    
    cur = conn.cursor()
    
    cur.execute("""
                SELECT      
                    x.start_node,
                    x.end_node,
                    ST_X(ST_AsText(ST_PointN(x.geom,1)))  start_x,
                    ST_Y(ST_AsText(ST_PointN(x.geom,1))) start_y,
                    ST_X(ST_AsText(ST_PointN(x.geom,2))) end_x,
                    ST_Y(ST_AsText(ST_PointN(x.geom,2))) end_y
                    FROM
                    (
                        SELECT     
                            a.start_node,
                            a.end_node,
                            a.geom
                        FROM  trips_topo.edge_data a  
                        --FROM  morris_test_20150114_new_topo.edge_data a  
                        GROUP BY a.edge_id,
                            a.start_node,
                            a.end_node,
                            a.geom
                    ) x    
                """)
    
    rows = cur.fetchall()
    cur.close()
    conn.close()
    G = g.multiGraph()

    
    for i in range(len(rows)):
        print i
        if i % 10000 == 0:
            print i
        
        a_coord_lat = Decimal(str(rows[i][2]))
        a_coord_lon = Decimal(str(rows[i][3]))
        b_coord_lat = Decimal(str(rows[i][4]))
        b_coord_lon = Decimal(str(rows[i][5]))
        
        #add edge to graph and pass on face information
        G.add_edge(rows[i][0],rows[i][1], (a_coord_lat, a_coord_lon), (b_coord_lat, b_coord_lon ))
    
    print 'DATA IMPORTED'
    G.n_maxID = G.increment_max_nodeID()
    return G