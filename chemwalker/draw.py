import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
import os, re
import networkx as nx
from chemwalker.visJS_module import visjs_network

def plotPannel(tlid, clusterid, score = 'chw_prob', nstruct=10):
    """Plots a pannel with candidate structures.
    Retrieves a list of candidate structures
    from each method and order by its score.
    Args:
        tabgnps: Node attributes table, with node order (do we need that?).
        tlid: a dict list with all node candidades (should we import the whole list?).
            to fetch.
        clusterid: cluster index.
        score: Scoring method to order the table.
        nstruct: Number of structures to display.
    Returns:
    A png pannel with the structures ordered by score, with structure id and score in
    the plot.
        example:
    Raises:
        IOError: ...
    """
    attr = tlid[tlid['cluster index']==clusterid]
    attr.sort_values(score, ascending=False, inplace=True)
    mols = [Chem.MolFromInchi(x) for x in attr["InChI"]]
    leg = []

    for i in attr.index[:nstruct]:
        cscore = str(round(float(attr.loc[i, score]), 3))
        identifier = attr.loc[i, "Identifier"]
        leg.append(f'{identifier}-{cscore}')

    img = Draw.MolsToGridImage(mols[:nstruct], molsPerRow=3,
                               subImgSize=(200,200),
                               legends=leg)
    return img

def drawSingleMol(inchi, fout):
    mol = Chem.MolFromInchi(str(inchi))
    Draw.MolToFile(mol, fout)

def plotGraph(dbmatch, tabgnps, tlid, net,
              method, dr, option=2, clusterid=None,
              comp=None, pos_scale=100):
    # option 1: connected component
    # option 2: direct neighbors

    if option==1:
        subnet = net[net['ComponentIndex']==comp]
        fname = f'comp_{comp}_{method}.html'
    elif option==2:
        c1 = net['CLUSTERID1']==clusterid
        c2 = net['CLUSTERID2']==clusterid
        subnet = net[c1 | c2]
        fname = f'clusterid_{clusterid}_{method}.html'

    G = nx.from_pandas_edgelist(subnet, source='CLUSTERID1',
                                 target='CLUSTERID2',
                                 edge_attr='Cosine')

    nodes = list(G.nodes())
    edges = list(G.edges())
    tabgnps.rename(columns={'cluster index': '#Scan#',
                            'parent mass': 'Precursor_MZ'},
                   inplace=True)

    ntabgnps = tabgnps.copy()[tabgnps['#Scan#'].isin(nodes)].reset_index(drop=True)
    dbmatch = pd.merge(ntabgnps[['#Scan#', 'Precursor_MZ']],
                       dbmatch[['#Scan#', 'Smiles', 'INCHI']],
                       on='#Scan#', how='left')
    dbmatch.fillna('', inplace=True)

    pos = nx.spring_layout(G)

    nodes_dict = []
    image = []
    for n in nodes:
        i = dbmatch[(dbmatch['#Scan#']==n) &
                    ((dbmatch.Smiles!='') |
                     (dbmatch.INCHI!='')) ].index.values
        img = ''
        if len(i):
            if dbmatch.loc[i[0], 'Smiles']!='' \
            and dbmatch.loc[i[0], 'INCHI']=='':
                smi = dbmatch.loc[i[0], 'Smiles']
                inchi = Chem.MolToInchi(Chem.MolFromSmiles(smi))
            else:
                inchi = dbmatch.loc[i[0], 'INCHI']
            img = 'LM_%s.png' % n
            fout = os.path.join(dr, img)

            pm = dbmatch.loc[i[0], 'Precursor_MZ']
            label = f'id:{n} m/z:{pm}'
        else:
            pm = dbmatch.loc[dbmatch['#Scan#']==n, 'Precursor_MZ'].values[0]
            label = f'id:{n} m/z:{pm}'

            tmp = tlid.copy()[tlid['cluster index']==n]
            if method=='MF' and len(tmp):
                inchi = tmp.iloc[0]['InChI']
                img = 'MF_%s.png' % n
            elif method=='RW' and len(tmp):
                tmp.sort_values(by='chw_prob',
                                ascending=False,
                                inplace=True)
                inchi = tmp.iloc[0]['InChI']
                img = 'RW_%s.png' % n

            fout = os.path.join(dr, img)

        try:
            drawSingleMol(inchi=inchi, fout=fout)
            image.append(img)
        except:
            image.append('')

        if img != '':
            shape = 'circularImage'
            address = fout
            ibool = True
            if 'MF' in address:
                nbcol = '#00FFFF'
            elif 'RW' in address:
                nbcol = '#0000FF'
            else:
                nbcol = '#00FF00'
            nbwidth = 20
        else:
            shape = 'triangle'
            address = 'undefined'
            ibool = False
            nbcol = '#FF0000'
            nbwidth = 1
        nodes_dict.append(
          {"id":n,
           #"degree":nx.degree(G,n),
           "degree": 1,
           "border_width": nbwidth,
           "label":label,
           "title":label,
           "node_size": 100,
           "node_shape_use_border_with_image":ibool,
           "border_color": nbcol,
           "node_shape": shape,
            "node_image":address,
            "x":pos[n][0]*pos_scale,
            "y":pos[n][1]*pos_scale
           }
        )

    node_map = dict(zip(nodes,range(len(nodes))))  # map to indices for source/target in edges

    # still can't change edge label
    elab = str(0.75)
    edges_dict = [{"source":node_map[edges[i][0]], "target":node_map[edges[i][1]],
                 "color":"gray","edge_label":0.75, "edge_width": 5} for i in range(len(edges))]

    f = open(fname, 'w+')
    f.write(
    # set some network-wide styles
    visjs_network(nodes_dict,edges_dict,
                  node_size_multiplier=50,
                  #node_color_border= '#0000FF',
                  #graph_width = '100%',
                  #graph_height = '100%',
                  node_size_transform = 'Math.abs',
                  node_font_size = 50,
                  edge_width=5,
                  navigation_buttons = True,
                  showButton = True,
                  output='html')
    )
    f.close()
