from chemwalker import rwalker
import pandas as pd
import json
import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

val = pd.read_csv('pcbi.1006089.s007.tsv', sep='\t')

task = '29d517e67067476bae97a32f2d4977e0'
net = pd.read_table('%s/net.tsv' % task)
tabgnps = pd.read_table('%s/tabgnps.tsv' % task)
with open('%s/lid.json' % task, 'r') as f:
    lid = json.load(f)

sval = pd.merge(tabgnps[['cluster.index', 'parent.mass']],
                val.loc[val['Processing batch']==1, ['cluster.index',
                                                     'parent.mass', 'MetFrag',
                                                     'Random', 'Fusion',
                                                     'Consensus']],
                left_on='cluster.index', right_on='cluster.index')


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

# position of metfrag
sval.head()
pd.DataFrame(lid[5]).head()
pd.DataFrame(lid[5])['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0]
np.where(pd.DataFrame(lid[5])['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0])

# see how changes of labels above match
os.listdir(task)
with open('%s/fusion.json' % task, 'r') as f:
    fusion = json.load(f)

np.where(pd.DataFrame(fusion[5]).sort_values(['fusion'], ascending=False)['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0])
sval.head()
with open('%s/consensus.json' % task, 'r') as f:
    consensus = json.load(f)
np.where(pd.DataFrame(consensus[5]).sort_values(['fusion2'], ascending=False)['InChIKey1']==val.loc[4, 'InChIKey Recovered'].split('-')[0])

fls = os.listdir('result_json')

seed = []
random = []
for fl in fls:
    with open('result_json/%s' % fl, 'r') as f:
        res = json.load(f)
    for k,v in json.loads(res).items():
        if 'seed' in k:
            seed.extend(v)
        else:
            random.extend(v)

seed = pd.DataFrame(seed)
random = pd.DataFrame(random)
rw = pd.merge(seed, random, left_on='cluster index', right_on='cluster index')
# why some are different?
sum(rw['metfrag_x']==rw['metfrag_y'])
rw= rw[rw['metfrag_x']==rw['metfrag_y']]

rw.columns = ['cluster index', 'metfrag_x', 'seed', 'metfrag_y', 'random']
sval = pd.merge(sval, rw[['cluster index', 'metfrag_x', 'seed']], left_on='cluster.index', right_on='cluster index')

sum(sval['MetFrag']==(sval['metfrag_x']+1))
sval = pd.merge(sval, rw[['cluster index', 'random']], left_on='cluster.index', right_on='cluster index')

sval['seed'] +=1
sval['random'] +=1

cumlist = {'MetFrag':[], 'Random':[], 'Fusion':[], 'Consensus':[], 'seed':[],
           'random':[]}
total = sval.shape[0]
for i in range(1, 21):
    cumlist['MetFrag'] = cumlist['MetFrag']+[sum(sval['MetFrag']<=i)/total]
    cumlist['Fusion'] = cumlist['Fusion']+[sum(sval['Random']<=i)/total]
    cumlist['Random'] = cumlist['Random']+[sum(sval['Consensus']<=i)/total]
    cumlist['Consensus'] = cumlist['Consensus']+[sum(sval['Fusion']<=i)/total]
    cumlist['seed'] = cumlist['seed']+[sum(sval['seed']<=i)/total]
    cumlist['random'] = cumlist['random']+[sum(sval['random']<=i)/total]


for k,v in cumlist.items():
    plt.plot(range(1, 21), v)


plt.legend(list(cumlist.keys()))
plt.show()
