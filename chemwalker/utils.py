import os
import re
import time
import json
import subprocess
import numpy as np
import pandas as pd
from chemwalker.rwalker import cand_pair, random_walk
from rdkit import Chem
import networkx as nx


precursor_ion_mode_positive = {
                             '[M+H]+': [1, 1.007276],
                             '[M+NH4]+': [18, 18.033823],
                             '[M+Na]+': [23, 22.989218],
                             '[M+K]+': [39, 38.963158],
                             '[M+CH3OH+H]+': [33, 33.033489],
                             '[M+ACN+H]+': [42, 42.033823],
                             '[M+ACN+Na]+': [64, 64.015765],
                             '[M+2ACN+H]+': [83, 83.060370],
                             '[M]+': [0, 0.0]
}
precursor_ion_mode_negative = {
                              '[M-H]-': [-1, -1.007276],
                              '[M+Cl]-': [35, 34.969402],
                              '[M+HCOO]-': [45, 44.998201],
                              '[M+CH3COO]-': [59, 59.013851],
                              '[M]-': [0, 0.0]
}

metfrag_param = {
                # data file containing mz intensity peak pairs (one per line)
                'PeakListPath': '',
                # database parameters -> how to retrieve candidates
                'MetFragDatabaseType': 'LocalPSV',
                'LocalDatabasePath': '',
                #NeutralPrecursorMolecularFormula = C9H11Cl3NO3PS
                #DatabaseSearchRelativeMassDeviation = PPM
                'NeutralPrecursorMass': 0,
                #IonizedPrecursorMass = 349.93356
                # peak matching parameters
                'FragmentPeakMatchAbsoluteMassDeviation': 0.01,
                'FragmentPeakMatchRelativeMassDeviation': 5,
                'PrecursorIonMode': 1,
                'IsPositiveIonMode': True,
                # scoring parameters
                'MetFragScoreTypes': 'FragmenterScore',
                'MetFragScoreWeights': 1.0,
                # output
                # SDF, XLS, CSV, ExtendedXLS, ExtendedFragmentsXLS
                'MetFragCandidateWriter': 'FragmentSmilesPSV',
                'SampleName': '',
                'ResultsPath': '.',
                # following parameteres can be kept as they are
                'MaximumTreeDepth': 2,
                'MetFragPreProcessingCandidateFilter': 'UnconnectedCompoundFilter',
                'MetFragPostProcessingCandidateFilter': 'InChIKeyFilter',
                'NumberThreads': 6
}

def filter_db(db, prmass, ppm, inchifilt=True):
    if inchifilt:
        db = db[~db.InChIKey1.duplicated()]
    cnames = ['MonoisotopicMass', 'InChI', 'Identifier', 'InChIKey2', 'InChIKey1', 'MolecularFormula']
    comp = 10**6 *((db.MonoisotopicMass-prmass).abs()/db.MonoisotopicMass)
    db = db[comp <= ppm].reset_index(drop=True)
    return db[cnames]

def run_metfrag(spectrum, db, cluster_index, adduct='[M+H]+', ppm=15, abs_diff=0.01,
                ispositive = True, metpath='MetFrag2.3-CL.jar'):
    if ispositive:
        prmass = spectrum['params']['pepmass'][0]-precursor_ion_mode_positive[adduct][1]
        metfrag_param['PrecursorIonMode'] = precursor_ion_mode_positive[adduct][0]
    else:
        metfrag_param['IsPositiveIonMode'] = False
        prmass = spectrum['params']['pepmass'][0]-precursor_ion_mode_negative[adduct][1]
        metfrag_param['PrecursorIonMode'] = precursor_ion_mode_negative[adduct][0]

    spec = zip(spectrum['m/z array'], spectrum['intensity array'])
    spec = ['%s\t%s\n' % x for x in spec]

    sname = 'spec_%s_data.txt' % spectrum['params']['scans']
    with open(sname, '+w') as f:
        for s in spec:
            f.write(s)

    if type(db)==str:
        metfrag_param['MetFragDatabaseType'] = db
    else:
        db_filt = filter_db(db, prmass, ppm)
        ndb = 'validation_db_%s.psv' % spectrum['params']['scans']
        db_filt.to_csv(ndb, sep='|', index=None)
        metfrag_param['LocalDatabasePath'] = ndb

    cname = 'cand_%s' % spectrum['params']['scans']
    metfrag_param['SampleName'] = cname
    metfrag_param['PeakListPath'] = sname
    metfrag_param['FragmentPeakMatchRelativeMassDeviation'] = ppm
    metfrag_param['FragmentPeakMatchAbsoluteMassDeviation'] = abs_diff
    metfrag_param['NeutralPrecursorMass'] = prmass

    pname = 'parameter_file_%s.txt' % spectrum['params']['scans']

    with open(pname, '+w') as f:
        for k,v in metfrag_param.items():
            f.write(f'{k} = {v}\n')

    subprocess.call(['java', '-jar', metpath, pname])
    cname = 'cand_%s.psv' % spectrum['params']['scans']
    try:
        metres = pd.read_csv(cname, sep='|')
        metres['cluster index'] = cluster_index
    except:
        metres = pd.DataFrame()

    os.remove(sname)
    if type(db)!=str:
        os.remove(ndb)
    os.remove(cname)
    os.remove(pname)

    return metres

def walk_conn_comp(net, spectra, tabgnps, dbmatch, comp_index, db,
                   method = 'RDKit7-linear', adduct='[M+H]+',
                   ppm=15, ispositive = True, metfrag_res=''):
    snet = net[net['ComponentIndex']==comp_index]

    # Obtain the nodes in the component
    nds = list(set(snet['CLUSTERID1'].tolist()+snet['CLUSTERID2'].tolist()))
    print('Component with %s nodes' % len(nds))

    # pd.DataFrame([[1, 1, 2, 2, 2], [0.9, 0.8, 0.75, 0.75, 0.6], ['A', 'A', 'B', 'B', 'B']], index=['cluster index', 'MQScore', 'PI']).T 
    dbmatch.rename(columns={'#Scan#': 'cluster index', 'Precursor_MZ': 'parent mass',
                            'RT_Query': 'RTMean', 'Compound_Name': 'LibraryID'}, inplace=True)
    mxscore = dbmatch.groupby('cluster index')['MQScore'].idxmax().values
    dbmatch = dbmatch.loc[mxscore]

    otabgnps = dbmatch.loc[dbmatch['cluster index'].isin(nds),
                       ['cluster index', 'parent mass', 'RTMean',
                        'LibraryID', 'Smiles', 'INCHI' ]]
    otabgnps.fillna('', inplace=True)
    idx_inchi = otabgnps[(otabgnps['Smiles']!='') &
                         (otabgnps['INCHI']=='')].index
    if len(idx_inchi):
        for i in idx_inchi:
            try:
                otabgnps.loc[i, 'INCHI'] = Chem.MolToInchi(Chem.MolFromSmiles(otabgnps.loc[i, 'Smiles']))
            except:
                pass
    inchikey = [Chem.InchiToInchiKey(x) if x!='' else '' for x in otabgnps['INCHI']]

    # Record the first block of InChIKey
    otabgnps['InChIKey1'] = [x.split('-')[0] if x!='' else '' for x in inchikey]

    lid = []

    assert len(spectra) == len(tabgnps),"Mismatch between spectra and attributes."

    start = time.time()
    print('Calculating in silico fragmentation with MetFrag...')
    for i in nds:
        if any((otabgnps['cluster index']==i) & (otabgnps['InChIKey1']!='')):
            continue
        j = np.where(tabgnps['cluster index']==i)[0][0]
        lid.append(run_metfrag(spectra[j], db, i,
                               adduct=adduct, ppm=ppm,
                               ispositive=ispositive))
    end = time.time()
    print('in silico fragmentation done in:', end - start, 'seconds')

    tlid = pd.concat(lid)
    if not len(tlid):
        raise ValueError("No candidate found by MetFrag")

    if  metfrag_res !='':
        with open('%s.json' % metrag_res) as f:
            json.dump(tlid.to_dict(), f)

    nr = (otabgnps['InChIKey1']!='').sum()
    if nr:
        cn = tlid.columns.tolist()
        dsource = pd.DataFrame([['']*len(cn)]*nr, columns=cn)
        dsource['InChI'] = otabgnps.loc[(otabgnps['InChIKey1']!=''), 'INCHI'].tolist()
        dsource['cluster index'] = otabgnps.loc[(otabgnps['InChIKey1']!=''), 'cluster index'].tolist()
        dsource['InChIKey1'] = otabgnps.loc[(otabgnps['InChIKey1']!=''), 'InChIKey1'].tolist()
        dsource['Identifier'] = otabgnps.loc[(otabgnps['InChIKey1']!=''), 'LibraryID'].tolist()
        dsource['Score'] = 1
        dsource['Score'] = dsource['Score'].astype(float)
        dsource['cluster index'] = dsource['cluster index'].astype(int)
        tlid = pd.concat([tlid, dsource]).reset_index(drop=True)

    start = time.time()
    print('Calculating pairwise candidate similarities...')
    scandpair = cand_pair(snet, tlid, method)
    end = time.time()
    print('similarities done in:', end - start, 'seconds')

    # include nodes without candidates?
    # and remove candidates for gnps ids
    miss_nds = set(tabgnps.loc[tabgnps['cluster index'].isin(nds), 'cluster index'])-set(tlid['cluster index'])
    if miss_nds:
        for n in miss_nds:
            npairs = snet.loc[(snet.CLUSTERID1==n) | (snet.CLUSTERID2==n),
                              ['CLUSTERID1', 'CLUSTERID2']].stack().unique().tolist()
            npairs = set(npairs)-{n}
            for p in npairs:
                pstr = '%s_' % p
                ploc = (scandpair.iloc[:,0].str.contains(pstr) |
                        scandpair.iloc[:,1].str.contains(pstr))
                # mean or min?
                if not ploc.sum():
                    pmean = scandpair.iloc[:,2].quantile(0.1)
                    dftmp = pd.DataFrame(['%s_nomatch' % n, '%s_nomatch' % p]).T
                else:
                    pmean = scandpair[ploc].iloc[:,2].mean()
                    pairs = scandpair.loc[ploc].iloc[:,[0,1]].stack()
                    pairs = pairs[pairs.str.contains(pstr)].unique().tolist()
                    dftmp = pd.DataFrame(zip(pairs, ['%s_nomatch' % n]*len(pairs)))
                dftmp[2] = pmean
                dftmp.columns = scandpair.columns
                scandpair = pd.concat([scandpair,
                                       dftmp]).reset_index(drop=True)

    G = nx.Graph()
    edge_list = scandpair.apply(lambda a: tuple(a), axis=1).tolist()
    G.add_weighted_edges_from(edge_list)

    glib = otabgnps.loc[otabgnps['InChIKey1']!='', 'cluster index'].index.tolist()
    source = []
    for i in glib:
        k,g = otabgnps.loc[i, ['InChIKey1', 'cluster index']].tolist()
        ksel = tlid[(tlid['cluster index']==g) & (tlid['InChIKey1']==k)]
        if len(ksel):
            idx = ksel['Identifier'].tolist()[0]
            print(f'Seed - InChIKey1:{k}, cluster index:{g}, Identifier:{idx}')
        else:
            print(f'No seed for - InChIKey1:{k}, cluster index:{g}')
        source.extend([x for x in G.nodes() if bool(re.search('%s_%s' % (g,idx), x))])

    if not len(source):
        raise ValueError("No GNPS id to propagate from")

    start = time.time()
    print('Walking on the graph...')
    p_t = random_walk(G, source)
    end = time.time()
    print('walking done in:', end - start, 'seconds')

    tlid['uid'] = tlid.apply(lambda a: f'{a["cluster index"]}_{a["Identifier"]}', axis=1)

    mol_probs = zip(G.nodes(), p_t.tolist())
    dprob = pd.DataFrame(list(mol_probs))
    dprob['cluster index'] = dprob.apply(lambda a: int(a[0].split('_')[0]), axis=1)
    dprob['Identifier'] = dprob.apply(lambda a: a[0].split('_')[1], axis=1)
    dprob.sort_values(['cluster index'], inplace=True)
    dprob.rename(columns={0: 'uid', 1: 'chw_prob'}, inplace=True)
    dprob.reset_index(inplace=True, drop=True)
    nprob = dprob.groupby('cluster index').apply(lambda grp: grp.chw_prob/grp.chw_prob.max())
    dprob['chw_prob'] = pd.DataFrame(nprob).reset_index()['chw_prob']
    tlid = pd.merge(tlid, dprob[['uid', 'chw_prob']], on='uid')
    return tlid

def val_known(tlid, dbmatch):
    dbmatch.rename(columns={'#Scan#': 'cluster index', 'Precursor_MZ': 'parent mass',
                            'RT_Query': 'RTMean', 'Compound_Name': 'LibraryID'}, inplace=True)
    mxscore = dbmatch.groupby('cluster index')['MQScore'].idxmax().values
    dbmatch = dbmatch.loc[mxscore]

    otabgnps = dbmatch.loc[dbmatch['cluster index'].isin(tlid['cluster index']),
                       ['cluster index', 'parent mass', 'RTMean',
                        'LibraryID', 'Smiles', 'INCHI' ]]
    otabgnps.fillna('', inplace=True)
    idx_inchi = otabgnps[(otabgnps['Smiles']!='') &
                         (otabgnps['INCHI']=='')].index
    if len(idx_inchi):
        for i in idx_inchi:
            try:
                otabgnps.loc[i, 'INCHI'] = Chem.MolToInchi(Chem.MolFromSmiles(otabgnps.loc[i, 'Smiles']))
            except:
                pass
    inchikey = [Chem.InchiToInchiKey(x) if x!='' else '' for x in otabgnps['INCHI']]

    # Record the first block of InChIKey
    otabgnps['InChIKey1'] = [x.split('-')[0] if x!='' else '' for x in inchikey]

    lrank = []
    for idx in otabgnps.index:
        if otabgnps.loc[idx, 'InChIKey1']!='':
            tmpid = tlid[tlid['cluster index']==otabgnps.loc[idx, 'cluster index']]
            if tmpid.shape[0]==0:
                continue
            mp = np.where(tmpid['InChIKey1']==otabgnps.loc[idx, 'InChIKey1'])[0]
            if len(mp)==0:
                continue

            tprob = tmpid.sort_values(['chw_prob'], ascending=False)
            rp = np.where(tprob['InChIKey1']==otabgnps.loc[idx, 'InChIKey1'])[0]
            if len(rp)==0:
                continue
            lrank.append({'cluster index': otabgnps.loc[idx, 'cluster index'], 'metfrag': mp[0], 'rw': rp[0]})
    return lrank


