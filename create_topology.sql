-- Create a topology
SELECT topology.CreateTopology('topology_name', 4326);
-- Drop a topology
SELECT topology.DropTopology('topology_name');
-- Add a topology column to your dataset containing the linestrings
SELECT topology.AddTopoGeometryColumn('topology_name', 'public', 'linestring_table_data', 'topo_geom', 'LINESTRING');
-- start the topology calculation, tolerance not needed, for more info check here:
-- http://blog.mathieu-leplatre.info/use-postgis-topologies-to-clean-up-road-networks.html
UPDATE linestring_table_data SET topo_geom = topology.toTopoGeom(geom, 'topology_name', 1);