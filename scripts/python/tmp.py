#!/usr/bin/env python

from openeye.oechem import *
import sys,os
import glob,psycopg2,traceback
import psycopg2.extras

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen2',user='medchem',host='javelin',password='medchem')

    def getMolFromBio(self, bio_number):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select m from rdk.mols where bio_number = %s",(bio_number,))
            result = cursor.fetchone()
            if result is not None:
                smiles = result['m']
                mol = OEGraphMol()
                ifs = oemolistream()
                ifs.SetFormat(OEFormat_SMI)
                ifs.openstring(smiles)
                OEReadMolecule(ifs,mol)
                return mol

            else:
                return None
        except:
            cursor.execute("rollback")
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()

files = glob.glob("*.txt")
db = BiogenDb()
for f in files:
    sdf_fname = "%s.sdf"%os.path.basename(f).split(".")[0]
    ofs = oemolostream()
    ofs.open(sdf_fname)
    print sdf_fname
    titles = []
    lines = open(f,"r").read().splitlines()
    for lineno,line in enumerate(lines):
        if lineno == 0:
            titles = line.split("\t")
            continue
        args = line.split("\t")
        mol = db.getMolFromBio(args[0])
        for id,arg in enumerate(args[1:]):
            OESetSDData(mol,titles[id+1],arg)
        OEWriteMolecule(ofs,mol)
    ofs.close()


