#!/usr/bin/env python

import sys,traceback,psycopg2,os
import psycopg2.extras
from openeye.oechem import *

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='frontier',user='medchem',host='javelin',password='medchem')

    def getASDIFromSmi(self, smiles):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select a.acdno from frontier_chemicals a inner join rdk.mols b on a.chemical_id=b.chemical_id where b.m@>%s and b.m<@%s",(smiles,smiles))
            result = cursor.fetchone()
            if result is not None:
                acdno = result['acdno']
                return acdno

            else:
                return None
        except:
            if cursor is not None:
                cursor.execute("rollback")
            print >>sys.stderr, smiles +"\thas error."
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__=="__main__":
    if len(sys.argv)!=3:
        print "Read in input.sdf, add asdi to output.sdf"
        print "Usage:%s input.sdf output.sdf"%os.path.basename(sys.argv[0])
    else:
        ifs  = oemolistream()
        ifs.open(sys.argv[1])
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        mol = OEGraphMol()
        db = BiogenDb()
        while OEReadMolecule(ifs,mol):
            smiles = OEMolToSmiles(mol)
            acdno = db.getASDIFromSmi(smiles)
            if acdno is None:
                OESetSDData(mol,"ASDI_ID","Unknown")
            else:
                OESetSDData(mol,"ASDI_ID",acdno)
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()

