import pandas as pd
import collections
import time
from chemwalker.rwalker import *
from chemwalker.gnps import Proteosafe
from rdkit import Chem
import json
import os
import requests
import numpy as np
import networkx as nx

methods = [(lambda x:Chem.RDKFingerprint(x,maxPath=4),'RDKit4'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=5),'RDKit5'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=6),'RDKit6'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=7),'RDKit7'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=4,branchedPaths=False),'RDKit4-linear'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=5,branchedPaths=False),'RDKit5-linear'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=6,branchedPaths=False),'RDKit6-linear'),
           (lambda x:Chem.RDKFingerprint(x,maxPath=7,branchedPaths=False),'RDKit7-linear'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,0),'MFP0'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,1),'MFP1'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,2),'MFP2'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,3),'MFP3'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,0,useFeatures=True),'FeatMFP0'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,1,useFeatures=True),'FeatMFP1'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,2,useFeatures=True),'FeatMFP2'),
           (lambda x:rdMolDescriptors.GetMorganFingerprint(x,3,useFeatures=True),'FeatMFP3'),
           (lambda x:rdMolDescriptors.GetHashedMorganFingerprint(x,0),'MFP0-bits'),
           (lambda x:rdMolDescriptors.GetHashedMorganFingerprint(x,1),'MFP1-bits'),
           (lambda x:rdMolDescriptors.GetHashedMorganFingerprint(x,2),'MFP2-bits'),
           (lambda x:rdMolDescriptors.GetHashedMorganFingerprint(x,3),'MFP3-bits'),
           (lambda x:rdMolDescriptors.GetAtomPairFingerprint(x),'AP'),
       (lambda x:rdMolDescriptors.GetTopologicalTorsionFingerprint(x),'TT'),
       (lambda x:rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(x),'AP-bits'),
               (lambda x:rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(x),'TT-bits'),
               (lambda x:rdMolDescriptors.GetMACCSKeysFingerprint(x),'MACCS'),
               (lambda x:pyAvalonTools.GetAvalonFP(x,512),'Avalon-512'),
               (lambda x:pyAvalonTools.GetAvalonFP(x,1024),'Avalon-1024'),
       ]


def val_known(G, p_t, otabgnps):
    mol_probs = zip(G.nodes(), p_t.tolist())
    dprob = pd.DataFrame(list(mol_probs))
    dprob['cluster index'] = dprob.apply(lambda a: int(a[0].split('_')[0]), axis=1)
    dprob['Identifier'] = dprob.apply(lambda a: a[0].split('_')[1], axis=1)
    dprob.sort_values(['cluster index'], inplace=True)
    lrank = []
    for idx in otabgnps.index:
        if otabgnps.loc[idx, 'InChIKey1']!='':
            tmpid = tlid[tlid['cluster.index']==otabgnps.loc[idx, 'cluster.index']]
            # understand why
            if tmpid.shape[0]==0:
                continue
            tmpid.reset_index(drop=False, inplace=True)
            tprob = dprob[dprob['cluster index']==otabgnps.loc[idx, 'cluster.index']]
            #tprob.sort_values([1], inplace=True, ascending=False)
            tprob = tprob.sort_values([1], ascending=False)
            mp = np.where(tmpid['InChIKey1']==otabgnps.loc[idx, 'InChIKey1'])[0]
            if len(mp)==0:
                continue
            rp = np.where(tprob['Identifier']==tmpid.loc[mp,'Identifier'].tolist()[0])[0]
            if len(rp)==0:
                continue
            lrank.append({'cluster index': otabgnps.loc[idx, 'cluster.index'], 'metfrag': mp[0], 'rw': rp[0]})
    return lrank

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

tasklist = ['29d517e67067476bae97a32f2d4977e0', 'd270e79876cb48deb6aabd52a4fc647e',
            'e2125577fe2646129becc248b96d42ba', '81e01fe178d3424686079903d908b536',
            'daa546b038604e5f83eaafb811bd0313', '61c8a0d01309408f8ecceb5b31dab1a8',
            '60fe9f77b3d04789997bf19aa1a0a828', '53f8494ff9e8423697eebf4e98d287f0',
            'c93a840100ec49bdbb3c12e5ed1e4790',
            '5fd60b02f8ab4274bf45fd5b715b5e0b']

directory = 'large_result'

if not os.path.exists(directory):
    os.makedirs(directory)


for task in tasklist:
    nap_result = Proteosafe(task, 'nap')
    nap_result.get_nap()
    net = nap_result.net
    tabgnps = nap_result.tabgnps
    lid = nap_result.lid
    nset = list(set(net['V7']))
    clist = []
    for ns in nset:
#        if task=='d270e79876cb48deb6aabd52a4fc647e':
#            ns += 70
#            if ns not in nset:
#                continue
        snet = net[net['V7']==ns]
        nds = list(set(snet['V1'].tolist()+snet['V2'].tolist()))
        otabgnps = tabgnps.loc[tabgnps['cluster.index'].isin(nds), ['cluster.index', 'parent.mass', 'RTMean', 'LibraryID', 'Smiles', 'INCHI' ]]
        #otabgnps['InChIKey1'] = otabgnps['INCHI'].apply(lambda a: Chem.InchiToInchiKey(a).split('-')[0] if type(a)==str else '')
        inchikey = [ Chem.InchiToInchiKey(x) if type(x)==str else ''for x in otabgnps['INCHI']]
        otabgnps['InChIKey1'] = [x.split('-')[0] if type(x)==str else '' for x in inchikey]
        start = time.time()
        # code still not in functions
        stabgnps = otabgnps.copy()
        seed_ctr=0.1
        seed_ctr = np.ceil(sum(otabgnps['InChIKey1']!='')*seed_ctr)
        seed_ctr = int(seed_ctr)
        exid = np.where(stabgnps['InChIKey1']!='')[0]
        exid =  stabgnps.index[list(exid[seed_ctr:])]
        stabgnps.loc[exid, 'InChIKey1'] = ''
        nidx = stabgnps.index
        tlid = []
        for x in nidx:
           if 'x' not in lid[x]:
               tmp = pd.DataFrame(lid[x])
               tmp['cluster.index'] = stabgnps['cluster.index'][x]
               tlid.append(tmp)
        for i in range(stabgnps.shape[0]):
            ind = stabgnps.index[i]
            if stabgnps.loc[ind, 'InChIKey1']!='' and \
               sum(tlid[i]['InChIKey1']==stabgnps.loc[ind, 'InChIKey1']):
                tlid[i] = tlid[i][tlid[i]['InChIKey1']==stabgnps.loc[ind, 'InChIKey1']]
        tlid = pd.concat(tlid)
        for func,method in methods:
            #method = 'MFP2-bits'
            scandpair = cand_pair(snet, tlid, method)
            G = nx.Graph()
            edge_list = scandpair.apply(lambda a: tuple(a), axis=1).tolist()
            G.add_weighted_edges_from(edge_list)
            glib = stabgnps.loc[stabgnps['InChIKey1']!='', 'cluster.index'].tolist()
            source = []
            for g in glib:
                source.extend([x for x in G.nodes() if bool(re.search('^%s_' % g, x))])
            p_t = random_walk(G, source)
            lr1 = val_known(G, p_t, otabgnps)
            #lr1 = walker(otabgnps, snet, lid, ncors=15)[0]
            end = time.time()
            print('Time for method %s for component %s with %s nodes and %s edges.' % (method, ns, len(nds), snet.shape[0] ))
            print(end - start)
            with open('%s/%s_%s_%s.json' % (directory, task, ns, method), 'w') as fp:
                json.dump(json.dumps({'%s_seed' % ns: lr1}, cls=MyEncoder), fp)
