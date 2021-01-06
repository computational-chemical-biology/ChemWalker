import os
import subprocess
import numpy as np
import pandas as pd


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
    os.remove(ndb)
    os.remove(cname)
    os.remove(pname)

    return metres
