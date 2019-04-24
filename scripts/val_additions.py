import os
os.chdir('../chemwalker/')
from rwalker import *
os.chdir('../ChemWalker')
import time
import os
import networkx as nx

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



task = '29d517e67067476bae97a32f2d4977e0'
net = pd.read_table('%s/net.tsv' % task)
tabgnps = pd.read_table('%s/tabgnps.tsv' % task)
with open('%s/lid.json' % task, 'r') as f:
    lid = json.load(f)

ns = 10
snet = net[net['V7']==ns]
nds = list(set(snet['V1'].tolist()+snet['V2'].tolist()))
otabgnps = tabgnps.loc[tabgnps['cluster.index'].isin(nds), ['cluster.index', 'parent.mass', 'RTMean', 'LibraryID', 'Smiles', 'INCHI' ]]
inchikey = [ Chem.InchiToInchiKey(x) if type(x)==str else ''for x in otabgnps['INCHI']]
otabgnps['InChIKey1'] = [x.split('-')[0] if type(x)==str else '' for x in inchikey]

start = time.time()
lr1 = walker(otabgnps, snet, lid, ncors=0)[0]
end = time.time()
print(end-start)

sum(pd.DataFrame(lr1)['metfrag']<pd.DataFrame(lr1)['rw'])
sum(pd.DataFrame(lr1)['metfrag']>pd.DataFrame(lr1)['rw'])

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
       nid = np.where(tmp['InChIKey1']==otabgnps.loc[x, 'InChIKey1'])[0][0]
       if tmp.shape[0]>5:
           lts = list(set([0, 1, 2, 3, nid]))
           tmp = tmp.iloc[lts]
       tmp['cluster.index'] = stabgnps['cluster.index'][x]
       tlid.append(tmp)

for i in range(stabgnps.shape[0]):
    ind = stabgnps.index[i]
    if stabgnps.loc[ind, 'InChIKey1']!='' and \
       sum(tlid[i]['InChIKey1']==stabgnps.loc[ind, 'InChIKey1']):
        tlid[i] = tlid[i][tlid[i]['InChIKey1']==stabgnps.loc[ind, 'InChIKey1']]


tlid = pd.concat(tlid)
method = 'MFP2-bits'

scandpair = cand_pair(snet, tlid, method)

G = nx.Graph()
edge_list = scandpair.apply(lambda a: tuple(a), axis=1).tolist()
G.add_weighted_edges_from(edge_list)
glib = stabgnps.loc[stabgnps['InChIKey1']!='', 'cluster.index'].tolist()
source = []
for g in glib:
    source.extend([x for x in G.nodes() if bool(re.search('^%s_' % g, x))])

p_t = random_walk(G, source)

lr2 = val_known(G, p_t, otabgnps)

sum(pd.DataFrame(lr1)['metfrag']<pd.DataFrame(lr1)['rw'])
sum(pd.DataFrame(lr2)['metfrag']<pd.DataFrame(lr2)['rw'])

sum(pd.DataFrame(lr1)['metfrag']>pd.DataFrame(lr1)['rw'])
sum(pd.DataFrame(lr2)['metfrag']>pd.DataFrame(lr2)['rw'])


