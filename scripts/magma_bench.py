
import pandas as pd
import collections 
import time 
from benchmark_auto import walker 
from rdkit import Chem
import json
import os
import requests
import numpy as np
from pyteomics import mgf 
from MAGMa import *
import sqlite3 

def getSingleSpectrum(mgffile, scanindex):
    spectrum = {} 
    f = open(mgffile) 
    lines = np.array(f.readlines())
    f.close()
    bg = np.where(lines=='BEGIN IONS\n') 
    ed = np.where(lines=='END IONS\n') 
    p = np.where(lines=='SCANS='+str(scanindex)+'\n') 
    bgp = bg[0][np.where((bg[0]-p[0])>0)[0][0]-1] 
    edp = ed[0][np.where((ed[0]-p[0])>0)[0][0]]+1 
    msms1 = []
    for line in lines[bgp:edp]:
        try:
            mz, rt = line.split(' ')  
            msms1.append((float(mz), float(rt)))  
        except:
            pass
    spectrum[scanindex] = pd.DataFrame(msms1, columns=['Mass', 'Intensity'])
    return(spectrum)

def byte2df(mres):
    lres = mres[0].decode("utf-8").split('\n') 
    lres = [re.sub(' \w+=', ' ', x).split(' ') for x in lres] 
    dfres = pd.DataFrame(lres) 
    dfres = dfres[dfres[0]!=''] 
    if dfres.shape[1]==6:
        dfres.columns = ['smiles', 'score', 'name', 'refscore', 'formula', 'mim']
        dfres['score'] = dfres['score'].astype(float)
        dfres['score'] = 1-(dfres['score']/dfres['score'].max())
    else:
        dfres.columns = ['smiles', 'name', 'refscore', 'formula', 'mim']
    return dfres

									        raise TypeError
# https://stackoverflow.com/questions/12734517/json-dumping-a-dict-throws-typeerror-keys-must-be-a-string
def stringify_keys(d):
    """Convert a dict's keys to strings if they are not."""
    for key in d.keys():
        # check inner dict
        if isinstance(d[key], dict):
            value = stringify_keys(d[key])
        else:
            value = d[key]
        if not isinstance(key, str):
            try:
                d[str(key)] = value
            except Exception:
                try:
                    d[repr(key)] = value
                except Exception:
                    raise
            # delete old key
            del d[key]
    return d

task = '29d517e67067476bae97a32f2d4977e0'
net = pd.read_table('%s/net.tsv' % task) 
tabgnps = pd.read_table('%s/tabgnps.tsv' % task) 
with open('%s/lid.json' % task, 'r') as f: 
    lid = json.load(f) 

#mgf = getSingleSpectrum('%s/allspectra.mgf' % task, tabgnps['cluster.index'][0])

lid2 = []
for ind in range(tabgnps.shape[0]):
    with mgf.read('%s/allspectra.mgf' % task) as reader: 
        for spectrum in reader: 
            if int(spectrum['params']['scans'])==tabgnps['cluster.index'][ind]:
                inta = spectrum['intensity array']
                if any(inta==0):
                    spectrum['intensity array'] = spectrum['intensity array'][inta!=0] 
                    spectrum['m/z array'] = spectrum['m/z array'][inta!=0] 
                mgf.write(iter([spectrum]), output='tmp.mgf')  
                break
    if type(lid[ind])==dict:
        lid2.append(pd.DataFrame())
        continue
    if pd.DataFrame(lid[ind]).shape[0]>0:     
        pd.DataFrame(lid[ind])[['SMILES', 'InChIKey1']].to_csv('tmp.smiles', sep='\t', header=None, index=None) 
        mres = runMAGMa('tmp.mgf', 'tmp.smiles', n='5') 
        lid2.append(byte2df(mres))
    else:
        lid2.append(pd.DataFrame())
    os.remove('tmp.mgf') 
    os.remove('tmp.smiles') 
    
mydict = [] 

for x in lid2:
    dtmp = stringify_keys(x.to_dict())
    for key,value in dtmp.items():
        dtmp[key] = stringify_keys(value)
    mydict.append(dtmp)

with open('%s/lid2.json' % task, 'w') as fout:
    json.dump(mydict, fout)


# docker run --rm -v $PWD:/data nlesc/magma read_ms_data --ms_data_format mgf glutathione.mgf results.db
# docker run --rm -v $PWD:/data nlesc/magma annotate -s hmdb results.db

conn = sqlite3.connect('results.db') 
c = conn.cursor() 
c.execute("SELECT name FROM sqlite_master WHERE type='table';") 
c.execute("SELECT * from 'fragments';")
c.fetchall() 

frags = pd.read_sql_query("SELECT * from 'fragments';", conn) 
run = pd.read_sql_query("SELECT * from 'run';", conn) 
mols = pd.read_sql_query("SELECT * from 'molecules';", conn) 
scans = pd.read_sql_query("SELECT * from 'scans';", conn) 
reactions = pd.read_sql_query("SELECT * from 'reactions';", conn)
peaks = pd.read_sql_query("SELECT * from 'peaks';", conn) 
