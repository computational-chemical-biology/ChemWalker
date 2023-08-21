from rdkit import DataStructs
from rdkit import rdBase
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import AllChem
from rdkit import Chem

import pandas as pd
import json
import collections
import random
import re
import xmltodict

import sys
import numpy as np
import networkx as nx
# sklearn just to normalize?
from sklearn.preprocessing import normalize
#from joblib import Parallel, delayed 
import multiprocessing
import scipy
import matplotlib.pyplot as plt

# https://github.com/TuftsBCB/Walker/blob/master/walker.py
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


# https://github.com/TuftsBCB/Walker/blob/master/walker.py
def normalize_cols(matrix):
    """ Normalize the columns of the adjacency matrix """
    return normalize(matrix, norm='l1', axis=0)

def getTanimoto(ms, method):
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
    method = [x[0] for x in methods if x[1]==method][0]
    m1, m2 = ms
    fp1 = method(Chem.MolFromInchi(m1))
    fp2 = method(Chem.MolFromInchi(m2))
    return DataStructs.TanimotoSimilarity(fp1,fp2)

# Needs expansion to encompass different fingerprints and combinations?? 
def pairSimilarity(inchipair, metric=DataStructs.TanimotoSimilarity,
                   fptype='circular'):
    """ Calculates fingerprint similarity """
    ms = [Chem.MolFromInchi(x) for x in inchipair]
    if fptype=='circular':
        fps = [AllChem.GetMorganFingerprint(x,2) for x in ms]
    else:
        fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    return metric(fps[0],fps[1])

def sigmoid(x, a=-9, b=0.6):
    """ Calculates sigmoid function as described in J. Mass Spectrom. 2013, 48, 291–298"""
    return 1/(1+np.exp(a*(x-b)))

def sc(metfrag, cosine, tn, alpha=0.3):
    """ Calculates scoring function as described in J. Mass Spectrom. 2013, 48, 291–298"""
    return alpha*metfrag+(1-alpha)*np.sum(sigmoid(cosine*tn))

def edge_weight(paramTuple):
    """ Calculates edge weight based on scoring presented above, or fing.
    simil. """
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

def edge_weight2(paramTuple):
    """ Calculates edge weight based on scoring presented above, or fing.
    simil. """
    cos, c1, inchi1, clust1, sc1, c2, inchi2, clust2, sc2, meansc, method = paramTuple
    # Choose between structural similarity and weighted score
    if meansc:
        metfrag = np.mean(list(map(float, [sc1, sc2])))
        try:
            #tn = pairSimilarity([inchi1, inchi2])
            tn = getTanimoto([inchi1, inchi2], method)
            paramTuple.extend([sc(metfrag, cos, tn)])
        except Exception as exc:
            print(exc)
            paramTuple.extend([metfrag])
    else:
        paramTuple.extend([getTanimoto([inchi1, inchi2], method)])
    return paramTuple

def cand_pair(snet, tlid, method, parallel = True, meansc = True, ncors=0):
    if parallel:
        dfparam = pd.merge(snet[['CLUSTERID1', 'CLUSTERID2', 'Cosine']],
                           tlid[['Identifier', 'InChI', 'cluster index', 'Score']],
                           left_on='CLUSTERID1', right_on='cluster index', how='left')
        dfparam = pd.merge(dfparam, tlid[['Identifier', 'InChI', 'cluster index', 'Score']],
                           left_on='CLUSTERID2', right_on='cluster index', how='left')
        dfnull = dfparam['cluster index_x'].isnull() | dfparam['cluster index_y'].isnull()
        if sum(dfnull):
            dfparam= dfparam[np.invert(dfnull)]
        dfparam['meansc'] = meansc
        dfparam.drop(['CLUSTERID1', 'CLUSTERID2'], inplace=True, axis=1)
        dfparam['cluster index_x'] = dfparam['cluster index_x'].astype(int)
        dfparam['cluster index_y'] = dfparam['cluster index_y'].astype(int)
        dfparam['method'] = method
        dfparam = dfparam.values.tolist()
        if(ncors):
            pool = multiprocessing.Pool(ncors)
        else:
            pool = multiprocessing.Pool()
        candpair = pool.map(edge_weight2, dfparam)
        candpair = pd.DataFrame(candpair)
        pool.close()
        candpair.reset_index(drop=True, inplace=True)
        # Keep track of cluster id and candidate id                                                      
        candpair[12] = candpair.apply(lambda a: '_'.join(map(str, a[[3, 1]])), axis=1 )
        candpair[13] = candpair.apply(lambda a: '_'.join(map(str, a[[7, 5]])), axis=1 )
        scandpair = candpair[[12, 13, 11]]
    return scandpair

def random_walk(graph, seed, restart_prob=0.8, step=0, epsilon=0.000001,
                niter=10000, sparce_matrix = True):
    p_0 = set_up_p0(seed, graph)
    thresh = 1
    p_t = np.copy(p_0)
    #sys.getsizeof(sog_not_normalized)
    # print small matrix to inspect col normalization
    # plot network with 3 nodes
    if sparce_matrix:
        w = nx.to_scipy_sparse_array(graph)
    else:
        w = nx.to_numpy_array(graph)
    w = normalize_cols(w)
    c = 0
    while (thresh > epsilon and c < niter):
        c += 1
        # first, calculate p^(t + 1) from p^(t)
        #p_t_1 = calculate_next_p(p_t, p_0)
        # p^{t+1} = (1-r)Wp^t+rp^0
        # np.dot matrix multiplication
        # Understand the role of re-start
        if sparce_matrix:
            walk = w.dot(p_t)
        else:
            walk = np.squeeze(np.asarray(np.dot(w, p_t)))
        no_restart = walk * (1 - restart_prob)
        restart = p_0 * restart_prob
        p_t_1 = np.add(no_restart, restart)
        if step:
            return p_t_1
        # calculate L1 norm of difference between p^(t + 1) and p^(t),
        # for checking the convergence condition
        thresh = np.linalg.norm(np.subtract(p_t_1, p_t), 1)
        # then, set p^(t) = p^(t + 1), and loop again if necessary
        # no deep copy necessary here, we're just renaming p
        p_t = p_t_1
    return p_t_1

def net_plt(graph, node_weight, layout_file='', edge_saling_factor=5,
            node_saling_factor=10000, fname=''):
    edw = [x[2]['weight'] for x in list(graph.edges(data=True))]
    if layout_file=='':
        pos = nx.spring_layout(graph)
    else:
        with open(layout_file) as fd:
            doc = xmltodict.parse(fd.read())
        #pos = {}
        for nd in doc['graph']['node']:
            #pos[float(nd['@label'])] = np.array([nd['graphics']['@x'], nd['graphics']['@y']])
            graph.node[float(nd['@label'])]['pos'] = (float(nd['graphics']['@x']),
            float(nd['graphics']['@y']))
        pos=nx.get_node_attributes(graph,'pos')
    nx.draw(graph, pos)
    nx.draw_networkx_edges(graph, pos, width=edge_saling_factor*edw)
    nx.draw_networkx_nodes(graph, pos, node_size=node_saling_factor*node_weight)
    # labels
    nx.draw_networkx_labels(graph, pos, font_size=12, font_family='sans-serif')
    if fname!='':
        nx.write_graphml_lxml(graph, fname)
    plt.axis('off')
    #plt.show()
    return [plt, graph]

