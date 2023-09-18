#!/usr/bin/env python
__author__ = 'jfeng1'
import glob
import os,sys
import psycopg2
import psycopg2.extras
import traceback
from progressbar import ProgressBar
import numpy
from openeye.oechem import *
from openeye.oegraphsim import *


class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen',user='medchem',host='javelin',password='medchem')

    def findExactMatch(self, smiles):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select m,bio_number from rdk.mols where m@>%s and m<@%s",(smiles,smiles))
            result = cursor.fetchone()
            if result is not None:
                smiles = result['m']
                bio_number = result['bio_number']
                mol = OEGraphMol()
                ifs = oemolistream()
                ifs.SetFormat(OEFormat_SMI)
                ifs.openstring(smiles)
                OEReadMolecule(ifs,mol)
                mol.SetTitle(str(bio_number))
                return mol

            else:
                return None
        except:
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()


    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    mol = OEGraphMol()
    ifs = oemolistream()
    ifs.open("darkmatter.sdf")
    molDb = OEMolDatabase(ifs)
    nrmols = molDb.GetMaxMolIdx()
    biogenDB = BiogenDb()
    ofs = oemolostream()
    ofs.open("biogen_darkmatter.sdf")
    ofs.SetFormat(OEFormat_SDF)
    pbar = ProgressBar(maxval=nrmols).start()
    for idx in range(0, nrmols):
        pbar.update(idx)
        mol = OEGraphMol()
        if molDb.GetMolecule(mol, idx):
            smiles = OEMolToSmiles(mol)
            mol = biogenDB.findExactMatch(smiles)
            if mol is not None:
                OEWriteMolecule(ofs,mol)
    ifs.close()
    ofs.close()

