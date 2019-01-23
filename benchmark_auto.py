from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit import Chem

import pandas as pd
import json
import collections 
import random
import re

import sys
import numpy as np
import networkx as nx
from sklearn.preprocessing import normalize
#from joblib import Parallel, delayed 
import multiprocessing  
import scipy

def set_up_p0(source, G):
    """ Set up and return the 0th probability vector. """
    p_0 = [0] * G.number_of_nodes()
    for source_id in source:
        try:
            # matrix columns are in the same order as nodes in original nx
            # graph, so we can get the index of the source node from the OG
            #source_index = G.nodes().index(source_id)
            nodes_list = list(G.nodes())
            source_index = nodes_list.index(source_id)
            p_0[source_index] = 1 / float(len(source))
        except ValueError:
            sys.exit("Source node {} is not in original graph. Source: {}. Exiting.".format(
                      source_id, source))
    return np.array(p_0)



def normalize_cols(matrix):
    """ Normalize the columns of the adjacency matrix """
    return normalize(matrix, norm='l1', axis=0)

def pairSimilarity(inchipair, metric=DataStructs.TanimotoSimilarity, fptype='circular'):
    ms = [Chem.MolFromInchi(x) for x in inchipair]
    if fptype=='circular':
        fps = [AllChem.GetMorganFingerprint(x,2) for x in ms]
    else:
        fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    return metric(fps[0],fps[1])

def sigmoid(x, a=-9, b=0.6):
    return 1/(1+np.exp(a*(x-b)))

def sc(metfrag, cosine, tn, alpha=0.3):
    return alpha*metfrag+(1-alpha)*np.sum(sigmoid(cosine*tn))

def edge_weight(paramTuple):
    cos, c1, inchi1, clust1, sc1, c2, inchi2, clust2, sc2, meansc = paramTuple
    # Choose between structural similarity and weighted score
    if meansc:
        metfrag = np.mean(list(map(float, [sc1, sc2]))) 
        try:
            tn = pairSimilarity([inchi1, inchi2]) 
            paramTuple.extend([sc(metfrag, cos, tn)])
        except Exception as exc:
            print(exc)
            paramTuple.extend([metfrag])
    else:   
        paramTuple.extend([pairSimilarity([inchi1, inchi2])])
    return paramTuple 
  
def walker(otabgnps, snet, lid, meansc = True, sparce_matrix=True, parallel=True, fix=True, rcandidate=False, seed_ctr=0.1, ncors=0):
    stabgnps = otabgnps.copy()
    # Think of resampling as using ids as seeds and not removing from list
    if seed_ctr:
        seed_ctr = np.ceil(sum(otabgnps['InChIKey1']!='')*seed_ctr) 
        seed_ctr = int(seed_ctr) 
        # add option to fixate the seed
        if fix:
            exid = np.where(stabgnps['InChIKey1']!='')[0]
            exid =  stabgnps.index[list(exid[seed_ctr:])]
            stabgnps.loc[exid, 'InChIKey1'] = '' 
        else:
            nseed = sum(stabgnps['InChIKey1']!='') - seed_ctr
            exid = random.sample(stabgnps.index.values.tolist(), nseed)
            exid.sort() 
            stabgnps.loc[exid, 'InChIKey1'] = '' 
       
    nidx = stabgnps.index  
    
    tlid = [] 
    for x in nidx: 
       if 'x' not in lid[x]: 
           tmp = pd.DataFrame(lid[x]) 
           tmp['cluster.index'] = stabgnps['cluster.index'][x] 
           tlid.append(tmp) 
    
    # select in the candidate list or replace GNPS id 
    if not rcandidate: 
        for i in range(stabgnps.shape[0]):
            ind = stabgnps.index[i]
            if stabgnps.loc[ind, 'InChIKey1']!='' and \
               sum(tlid[i]['InChIKey1']==stabgnps.loc[ind, 'InChIKey1']):
                tlid[i] = tlid[i][tlid[i]['InChIKey1']==stabgnps.loc[ind, 'InChIKey1']] 
    
    tlid = pd.concat(tlid) 
    
    if parallel:    
        dfparam = pd.merge(snet[['V1', 'V2', 'V5']], tlid[['Identifier', 'InChI', 'cluster.index', 'Score']], left_on='V1', right_on='cluster.index', how='left')
        dfparam = pd.merge(dfparam, tlid[['Identifier', 'InChI', 'cluster.index', 'Score']], left_on='V2', right_on='cluster.index', how='left')
        dfnull = dfparam['cluster.index_x'].isnull() | dfparam['cluster.index_y'].isnull()
        if sum(dfnull):
            dfparam= dfparam[np.invert(dfnull)] 
        dfparam['meansc'] = meansc
        dfparam.drop(['V1', 'V2'], inplace=True, axis=1) 
        dfparam['cluster.index_x'] = dfparam['cluster.index_x'].astype(int) 
        dfparam['cluster.index_y'] = dfparam['cluster.index_y'].astype(int) 
        dfparam = dfparam.values.tolist() 
        
        if(ncors):
            pool = multiprocessing.Pool(ncors) 
        else:
            pool = multiprocessing.Pool() 
        candpair = pool.map(edge_weight, dfparam) 
        candpair = pd.DataFrame(candpair)
        candpair.reset_index(drop=True, inplace=True) 
        
        # Keep track of cluster id and candidate id 
        candpair[11] = candpair.apply(lambda a: '_'.join(map(str, a[[3, 1]])), axis=1 )
        candpair[12] = candpair.apply(lambda a: '_'.join(map(str, a[[7, 5]])), axis=1 )
        scandpair = candpair[[11, 12, 10]]
    else:
        #res = Parallel(n_jobs=2)(delayed(edge_weight)(x) for x in dfparam)
        # create all candidates edge list
        # can be done in parallel        
        candpair = []
        
        for idx in snet.index:
            cand1 = tlid[tlid['cluster.index']==snet.loc[idx, 'V1']][['cluster.index', 'InChI', 'Identifier', 'Score']].reset_index(drop=True)
            cand2 = tlid[tlid['cluster.index']==snet.loc[idx, 'V2']][['cluster.index', 'InChI', 'Identifier', 'Score']].reset_index(drop=True)
            candtmp = []
            for c1 in cand1.index:
                for c2 in cand2.index:
                     tmp = cand1.iloc[c1].tolist()+cand2.iloc[c2].tolist()  
                     # Choose between structural similarity and weighted score
                     if meansc:
                         metfrag = np.mean(list(map(float, [cand1.loc[c1, 'Score'], cand2.loc[c2, 'Score']]))) 
                         cos =  snet.loc[idx, 'V5']
                         tn = pairSimilarity([cand1.loc[c1, 'InChI'], cand2.loc[c2, 'InChI']]) 
                         tmp.extend([sc(metfrag, cos, tn)])
                     else:   
                         tmp.extend([pairSimilarity([cand1.loc[c1, 'InChI'], cand2.loc[c2, 'InChI']])])
                     candtmp.append(tmp)
            # normalize the scores by node pair?
            tcand = pd.DataFrame(candtmp)
            candpair.append(tcand)
        
        candpair = pd.concat(candpair)
        candpair.reset_index(drop=True, inplace=True) 
        
        # Keep track of cluster id and candidate id 
        candpair[9] = candpair.apply(lambda a: '_'.join(map(str, a[[0, 2]])), axis=1 )
        candpair[10] = candpair.apply(lambda a: '_'.join(map(str, a[[4, 6]])), axis=1 )
        scandpair = candpair[[9, 10, 8]]
           
    # Generate weigheted graph  
    G = nx.Graph()
    edge_list = scandpair.apply(lambda a: tuple(a), axis=1).tolist()
    G.add_weighted_edges_from(edge_list)
    
    # convergence criterion - when vector L1 norm drops below 10^(-6)
    # (this is the same as the original RWR paper)
    CONV_THRESHOLD = 0.000001
    
    #def run_exp(source, restart_prob, og_prob, node_list=[]):
    if sum(stabgnps['InChIKey1']!='')==0 or rcandidate:
        #ndegree = [x for _,x in sorted(zip(G.degree().values(),G.degree().keys()))]
        ndegree = [x for _,x in sorted(zip(dict(G.degree()).values(),dict(G.degree()).keys()))] 
        seed_ctr = np.ceil(len(ndegree)*seed_ctr) 
        source = ndegree[-1:]
    else:
        glib = stabgnps.loc[stabgnps['InChIKey1']!='', 'cluster.index'].tolist() 
        # Sample nodes when there is no library id
        source = []
        for g in glib:
            source.extend([x for x in G.nodes() if bool(re.search('^%s_' % g, x))])
    
    #source = ['0', '1']
    # 'Restart probability for random walk', default=0.7
    # Lower the re-start probability when it is random?
    restart_prob = 0.7 
    # '--original_graph_prob', default=0.1
    og_prob = 0.1
    # 'seed', help='Seed file, to pull start nodes from'
    node_list = []
    """ Run a multi-graph random walk experiment, and print results.
    Parameters:
    -----------
        source (list):        The source node indices (i.e. a list of Entrez
                              gene IDs)
        restart_prob (float): As above
        og_prob (float):      As above
    """
    restart_prob = restart_prob
    og_prob = og_prob
    # set up the starting probability vector
    #p_0 = np.array([0.5, 0.5, 0.0, 0.0, 0.0])
    p_0 = set_up_p0(source, G)
    diff_norm = 1
    # this needs to be a deep copy, since we're reusing p_0 later
    p_t = np.copy(p_0)
    
    """ Build column-normalized adjacency matrix for each graph.
    NOTE: these are column-normalized adjacency matrices (not nx
          graphs), used to compute each p-vector
    """
    #sys.getsizeof(sog_not_normalized)
    # print small matrix to inspect col normalization
    # plot network with 3 nodes
    if sparce_matrix:
        og_not_normalized = nx.to_scipy_sparse_matrix(G) 
    else:
        og_not_normalized = nx.to_numpy_matrix(G)
    
    og_matrix = normalize_cols(og_not_normalized)
    
    while (diff_norm > CONV_THRESHOLD):
        # first, calculate p^(t + 1) from p^(t)
        #p_t_1 = calculate_next_p(p_t, p_0)
        # p^{t+1} = (1-r)Wp^t+rp^0
        # np.dot matrix multiplication
        # Understand the role of re-start
        if sparce_matrix:
            epsilon = og_matrix.dot(p_t)
        else:
            epsilon = np.squeeze(np.asarray(np.dot(og_matrix, p_t)))
        no_restart = epsilon * (1 - restart_prob)
        restart = p_0 * restart_prob
        p_t_1 = np.add(no_restart, restart)
        # calculate L1 norm of difference between p^(t + 1) and p^(t),
        # for checking the convergence condition
        diff_norm = np.linalg.norm(np.subtract(p_t_1, p_t), 1)
        # then, set p^(t) = p^(t + 1), and loop again if necessary
        # no deep copy necessary here, we're just renaming p
        p_t = p_t_1
       
    # create probability for each node 
    # Name the data frame columns
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
       
    return [lrank, dprob]  

