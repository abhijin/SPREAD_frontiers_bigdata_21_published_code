DESC='''
Convert simulations data to CART
'''

import numpy as np
from os import environ,getenv
from os.path import basename
import pandas as pd
from pdb import set_trace
import sqlite3

# read table
conn = sqlite3.connect(getenv('RESULTS_DB'))
df = pd.read_sql('SELECT * FROM summary',conn).rename(columns = {
    'moore_range': 'range',
    'random_seed': 'seed'})
print(df.shape)
net_prop = pd.read_sql('SELECT * FROM mpsn_props',conn)
set_trace()
net_prop = net_prop[(net_prop.number_of_nodes == 4096) &
        (net_prop.number_of_regions == 16) &
        (net_prop.locality_graph == 'clique')][['range',
            'locality_size',
            'epsilon',
            'seed',
            'spectral_radius',
            'diameter']]
net_prop = net_prop.set_index(['range', 'locality_size', 'epsilon', 'seed'])
df['locality_size'] = df.network.str.replace('.*s_','').str.replace('/.*','').astype(int)
df['epsilon'] = df.network.str.replace('.*p_','').str.replace('/.*','').astype(float)
df = df.set_index(['range', 'locality_size', 'epsilon', 'seed'])

df = df.join(net_prop)

df[['alpha_S', 'alpha_L', 'alpha_LD', 'time_step', 
       'infections_mean', 'spectral_radius', 
       'diameter']].to_csv('simulation_data.csv')
