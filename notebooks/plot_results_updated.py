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
import matplotlib.pyplot as plt

val = pd.read_csv('data/pcbi.1006089.s007.tsv', sep='\t', skiprows=2)
val.head()

df_list = []
fls = os.listdir('large_result')
for fl in fls:
    with open(os.path.join('large_result/', fl)) as f:
        try:
            tmp = json.loads(json.load(f))
        except:
            continue
        tmp = pd.DataFrame(tmp[list(tmp.keys())[0]])
        _, n, fg = fl.split('_')
        tmp['comp'] = n
        tmp['fingerprint'] = fg.split('.')[0]
        df_list.append(tmp)

task = '29d517e67067476bae97a32f2d4977e0'
nap_result = Proteosafe(task, 'nap')
nap_result.get_nap()
net = nap_result.net
tabgnps = nap_result.tabgnps
sval = pd.merge(tabgnps[['cluster.index', 'parent.mass']],
                val.loc[val['Processing batch']==1, ['cluster.index',
                                                     'parent.mass', 'MetFrag',
                                                     'Random', 'Fusion',
                                                     'Consensus']],
                left_on='cluster.index', right_on='cluster.index')
sval.shape
cumlist = {'MetFrag':[], 'Random':[], 'Fusion':[], 'Consensus':[]}
total = sval.shape[0]
# labels got mixed during export of final table
for i in range(1, 21):
    cumlist['MetFrag'] = cumlist['MetFrag']+[sum(sval['MetFrag']<=i)/total]
    cumlist['Fusion'] = cumlist['Fusion']+[sum(sval['Random']<=i)/total]
    cumlist['Random'] = cumlist['Random']+[sum(sval['Consensus']<=i)/total]
    cumlist['Consensus'] = cumlist['Consensus']+[sum(sval['Fusion']<=i)/total]

for k,v in cumlist.items():
    plt.plot(range(1, 21), v)

plt.legend(list(cumlist.keys()),
           shadow=True, handlelength=1.5, fontsize=16)
plt.show()

lid = nap_result.lid
pd.DataFrame(lid[5]).head()
pd.DataFrame(lid[5])['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0]
np.where(pd.DataFrame(lid[5])['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0])
fusion = nap_result.fusion
np.where(pd.DataFrame(fusion[5]).sort_values(['fusion'], ascending=False)['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0])
tabgnps.head()

fgs = df.fingerprint.unique()
fdf = []
for f in fgs:
    sval = pd.merge(df[df.fingerprint==f],
                    val.loc[val['Processing batch']==1, ['cluster.index',
                                                         'parent.mass', 'MetFrag',
                                                         'Random', 'Fusion',
                                                         'Consensus']],
                    left_on='cluster index', right_on='cluster.index')
    sval['metfrag'] += 1
    sval['rw'] += 1
    sval = sval[sval['metfrag']==sval['MetFrag']]
    print(sval.shape[0], 'for fingerprint', f)
    cumlist = {'MetFrag':[], 'Random':[], 'Fusion':[], 'Consensus':[], 'seed':[]}
    total = sval.shape[0]
    for i in range(1, 21):
        cumlist['MetFrag'] = cumlist['MetFrag']+[sum(sval['MetFrag']<=i)/total]
        cumlist['Fusion'] = cumlist['Fusion']+[sum(sval['Random']<=i)/total]
        cumlist['Random'] = cumlist['Random']+[sum(sval['Consensus']<=i)/total]
        cumlist['Consensus'] = cumlist['Consensus']+[sum(sval['Fusion']<=i)/total]
        cumlist['seed'] = cumlist['seed']+[sum(sval['rw']<=i)/total]
    tmp = pd.DataFrame(cumlist)
    tmp['fingerprint'] = f
    fdf.append(tmp)


for f in fgs:
    plt.plot(range(1, 21), fdf.loc[fdf.fingerprint==f, 'seed'])


plt.plot(range(1, 21), fdf.loc[fdf.fingerprint==f, 'MetFrag'], linestyle='--')
plt.plot(range(1, 21), fdf.loc[fdf.fingerprint==f, 'Random'], linestyle='--')

plt.legend(fgs.tolist()+['MetFrag', 'Random'])
plt.show()
