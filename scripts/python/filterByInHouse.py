#!/usr/bin/env python

import sys,traceback,psycopg2,os
import psycopg2.extras
from openeye.oechem import *

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen2',user='medchem',host='javelin',password='medchem')

    def getBioNumberFromSmi(self, smiles):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select bio_number from rdk.mols where m@>%s and m<@%s",(smiles,smiles))
            result = cursor.fetchone()
            if result is not None:
                bio_number = result['bio_number']
                return bio_number

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
        print "Read in input.sdf, export all molecules that do not have a BIO_NUMBER to output.sdf"
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
            bio_number = db.getBioNumberFromSmi(smiles)
            if bio_number is None:
                OEWriteMolecule(ofs,mol)
            else:
                print >>sys.stderr,"%s %s"%(smiles,bio_number)
        ifs.close()
        ofs.close()

