#!/usr/bin/env python

import sys,traceback,psycopg2,os
import psycopg2.extras
from openeye.oechem import *

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='aldrich',user='medchem',host='javelin',password='medchem')

    def getAldrichIdFromSmi(self, smiles):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select structure_id from rdk.mols where m@>%s and m<@%s",(smiles,smiles))
            result = cursor.fetchone()
            if result is not None:
                structure_id = result['structure_id']
                return structure_id

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
        print "Read in input.sdf, export all molecules that have aldrich id to output.sdf"
        print "Usage:%s input.sdf output.sdf"%os.path.basename(sys.argv[0])
    else:
        ifs  = oemolistream()
        ifs.open(sys.argv[1])
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        mol = OEGraphMol()
        db = BiogenDb()
        print >> sys.stderr,"Found:"
        while OEReadMolecule(ifs,mol):
            smiles = OEMolToSmiles(mol)
            structure_id = db.getAldrichIdFromSmi(smiles)
            if structure_id is not None:
                OESetSDData(mol,"AldrichId",structure_id)
                OEWriteMolecule(ofs,mol)
                print >>sys.stderr,"%s %s"%(smiles,structure_id)
        ifs.close()
        ofs.close()

