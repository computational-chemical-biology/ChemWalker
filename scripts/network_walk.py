#!/home/ridasilva/miniconda3/envs/eic/bin/python

import pandas as pd
import collections
import time
from benchmark_auto import walker
from rdkit import Chem
import json
import os
import requests
import numpy as np

class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(MyEncoder, self).default(obj)

tasklist = ['29d517e67067476bae97a32f2d4977e0', 'd270e79876cb48deb6aabd52a4fc647e', 'e2125577fe2646129becc248b96d42ba', '81e01fe178d3424686079903d908b536', 'daa546b038604e5f83eaafb811bd0313', '61c8a0d01309408f8ecceb5b31dab1a8', '60fe9f77b3d04789997bf19aa1a0a828', '53f8494ff9e8423697eebf4e98d287f0', 'c93a840100ec49bdbb3c12e5ed1e4790', '5fd60b02f8ab4274bf45fd5b715b5e0b']

os.chdir('/home/ridasilva/walker/full_benchmark')

for task in tasklist[1:]:
    net = pd.read_table('%s/net.tsv' % task)
    tabgnps = pd.read_table('%s/tabgnps.tsv' % task)
    with open('%s/lid.json' % task, 'r') as f:
        lid = json.load(f)

    nset = list(set(net['V7']))

    clist = []
    for ns in nset:
        if task=='d270e79876cb48deb6aabd52a4fc647e':
            ns += 70
            if ns not in nset:
                continue
        snet = net[net['V7']==ns]
        nds = list(set(snet['V1'].tolist()+snet['V2'].tolist()))
        otabgnps = tabgnps.loc[tabgnps['cluster.index'].isin(nds), ['cluster.index', 'parent.mass', 'RTMean', 'LibraryID', 'Smiles', 'INCHI' ]]
        #otabgnps['InChIKey1'] = otabgnps['INCHI'].apply(lambda a: Chem.InchiToInchiKey(a).split('-')[0] if type(a)==str else '')
        inchikey = [ Chem.InchiToInchiKey(x) if type(x)==str else ''for x in otabgnps['INCHI']]
        otabgnps['InChIKey1'] = [x.split('-')[0] if type(x)==str else '' for x in inchikey]
        start = time.time()
        lr1 = walker(otabgnps, snet, lid, ncors=15)[0]
        end = time.time()
        print('Time for component %s with %s nodes and %s edges.' % (ns, len(nds), snet.shape[0] ))
        print(end - start)
        with open('%s_%s_%s.json' % (task, ns, 'seed'), 'w') as fp:
            json.dump(json.dumps({'%s_seed' % ns: lr1}, cls=MyEncoder), fp)
        start = time.time()
        lr2 = walker(otabgnps, snet, lid, rcandidate=True, ncors=15)[0]
        end = time.time()
        print('Time for component %s with %s nodes and %s edges with random nodes.' % (ns, len(nds), snet.shape[0] ))
        print(end - start)
        with open('%s_%s_%s.json' % (task, ns, 'random'), 'w') as fp:
            json.dump(json.dumps({'%s_random' % ns: lr2}, cls=MyEncoder), fp)

