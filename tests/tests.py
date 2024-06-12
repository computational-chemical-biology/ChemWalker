import unittest
import os
from chemwalker.gnps import Proteosafe
from chemwalker.utils import *
import pandas as pd
import numpy as np

class TestWalk(unittest.TestCase):
    def setUp(self):
        taskid = 'b6d12e1e42b64ad6b20e7d38d5a4214b'
        gnps_result = Proteosafe(taskid, 'FBMN-gnps2')
        gnps_result.get_gnps()
        self.gnps_result = gnps_result
        db = get_db()
        self.db = db

    def test_gnps_api(self):
        n = self.gnps_result.net.shape[0]
        self.assertEqual(n, 580, f'Edge list with {n} rather than 580 edges')

    def test_comp(self):
        nds = self.gnps_result.gnps.component==11
        pos = self.gnps_result.dbmatch['#Scan#'].isin(self.gnps_result.gnps.loc[nds, 'cluster index'])
        #n = (self.gnps_result.dbmatch.loc[pos, 'Smiles'].notnull() | self.gnps_result.dbmatch.loc[pos, 'INCHI'].notnull()).sum()
        self.assertEqual(nds.sum(), 8, f'Number of nodes {nds.sum()} in component 11 rather than 8 edges')

    def test_db(self):
        n = self.db.shape[0]
        self.assertEqual(n, 367204, f'Number of entries {n} in db rather than 367204 compounds')

    def test_metfrag(self):
        metpath = os.path.abspath('bin/MetFrag2.3-CL.jar')
        res = run_metfrag(spectrum=self.gnps_result.spectra[0],
                          db=self.db, cluster_index=5, **{'metpath': metpath})
        n = res.shape[0]
        self.assertEqual(n, 54, f'Number of entries {n} in MetFrag results rather than 54 compounds')

    def test_random_walk(self):
        metpath = os.path.abspath('bin/MetFrag2.3-CL.jar')
        tlid = walk_conn_comp(net=self.gnps_result.net, spectra=self.gnps_result.spectra, tabgnps=self.gnps_result.gnps,
                              dbmatch=self.gnps_result.dbmatch.copy(), db=self.db, comp_index=11,
                              metpath=metpath)
        n = tlid.shape[0]
        self.assertEqual(n, 67, f'Number of entries {n} in RW results rather than 67 compounds')



if __name__ == '__main__':
    unittest.main()

