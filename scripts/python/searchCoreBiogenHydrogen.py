#!/usr/bin/env python
from openeye.oechem import *
from openeye.oemedchem import *
import sys, os
import psycopg2
import psycopg2.extras
import traceback
class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen2',user='medchem',host='javelin.corp.biogen.com',password='medchem')


    def getSmilesBySubstructure(self,sub_smi):
        smiList = []
        cursor = None
        try:
            print "finding "+sub_smi
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select m,bio_number from rdk.mols where m@>%s",(sub_smi,))
            results = cursor.fetchall()
            if results is not None:
                for r in results:
                    smiles = r['m']
                    bio_number = r['bio_number']
                    smiList.append("%s %s"%(smiles,bio_number))
                return smiList
            else:
                return None
        except:
            if cursor is not None:
                cursor.execute("rollback")
            print >>sys.stderr, sub_smi +"\thas error."
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def getSmilesBySmarts(self,sub_smi):
        smiList = []
        cursor = None
        try:
            sub_smi = sub_smi.replace("*","[*]").replace("[nH]","n")
            print "finding "+sub_smi
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select m,bio_number from rdk.mols where m@>%s::qmol",(sub_smi,))
            results = cursor.fetchall()
            if results is not None:
                for r in results:
                    smiles = r['m']
                    bio_number = r['bio_number']
                    smiList.append("%s %s"%(smiles,bio_number))
                return smiList
            else:
                return None
        except:
            if cursor is not None:
                cursor.execute("rollback")
            print >>sys.stderr, sub_smi +"\thas error."
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

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


if __name__ == "__main__":
    if len(sys.argv)!=3:
        print "Convert compounds to core, and search over BiogenDB"
        print "%s input.sdf output.sdf"%sys.argv[0]
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        dict = {}
        db = BiogenDb()
        coreKeys = []
        while OEReadMolecule(ifs,mol):
            newmol = OEGraphMol(mol)
            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol)
            OEAddExplicitHydrogens(newmol)
            smiles = OECreateSmiString(newmol,OESMILESFlag_ImpHCount)
            print(smiles)
            smiList = db.getSmilesBySmarts(smiles)
            OESetSDData(newmol,"parentSmiles",smiles)
            OESetSDData(newmol,"parentHits","%d"%len(smiList))

            OEWriteMolecule(ofs,newmol)

        ifs.close()
        ofs.close()
