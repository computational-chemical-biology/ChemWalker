from chemwalker.gnps import Proteosafe
from chemwalker.utils import run_metfrag, walk_conn_comp, val_known
from chemwalker.draw import *
import pandas as pd
import numpy as np

import click

@click.group()
def chemwalker():
    pass

@chemwalker.command()
@click.option("--taskid",
              help="GNPS task id")
@click.option("--workflow",
              default='FBMN',
              help="Workflow type, either FBMN or V2")
@click.option("--comp",
              help="Component (Molecular family) index")
@click.option("--db",
              help="Database file")
@click.option("--out",
              help="Output file name")
def random_walk(taskid, workflow, comp, db, out):
    gnps_result = Proteosafe(taskid, workflow)
    gnps_result.get_gnps()
    net = gnps_result.net
    gnps_tab = gnps_result.gnps
    spectra = gnps_result.spectra
    match_tab = gnps_result.dbmatch

    db = pd.read_csv(db, sep='|')

    tlid = walk_conn_comp(net=net, spectra=spectra, tabgnps=gnps_tab,
                          dbmatch=match_tab, db=db, comp_index=comp)

    tlid.to_csv('%s.tsv' % out, sep='\t', index=None)

if __name__ == '__main__':
    chemwalker()
