#!/bin/bash
sqlite3 ../results/results.db << 'EOF'
CREATE TABLE IF NOT EXISTS mpsn_props (
number_of_nodes INTEGER,
number_of_regions INTEGER,
range REAL,
locality_size INTEGER,
locality_graph TEXT,
long_distance_type TEXT,
epsilon REAL,
seed INTEGER,
number_of_edges INTEGER,
GS_number_of_edges INTEGER,
GL_number_of_edges INTEGER,
GLD_number_of_edges INTEGER,
max_core_number INTEGER,
spectral_radius REAL,
laplacian REAL,
diameter INTEGER,
PRIMARY KEY(number_of_nodes,number_of_regions,range,locality_size,locality_graph,
long_distance_type,epsilon,seed));
EOF
