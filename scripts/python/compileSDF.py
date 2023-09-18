#!/usr/bin/env python

import glob
import os,sys
import psycopg2
import psycopg2.extras
import traceback
import numpy
from openeye.oechem import *
from openeye.oegraphsim import *
DIRECTORY = "/Users/jfeng1/Datasets/MOA/"

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
            cursor.execute("rollback")
            print smiles +"\thas error."
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    db = BiogenDb()
    # print db.getBioNumberFromSmi("n2(c1nc(nc(c1nc2)NCc3ccccc3)NC(CO)CC)C(C)C")
    molDict = {}
    vendorDict = {}
    bioNumberDict = {}
    path = os.path.join(DIRECTORY,"*.sdf")
    files = glob.glob(path)
    for f in files:
        fname = os.path.basename(f)
        vendor_name = fname.split(".")[0]
        ifs = oemolistream()
        ifs.open(fname)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            mol1 = OEGraphMol(mol)
            OEDeleteEverythingExceptTheFirstLargestComponent(mol)
            smiles = OEMolToSmiles(mol)
            if not bioNumberDict.has_key(smiles):
                bio_number = db.getBioNumberFromSmi(smiles)
                if bio_number is not None:
                    bioNumberDict[smiles] = bio_number

            if not molDict.has_key(smiles):
                molDict[smiles] = mol1
                vendorDict[smiles] = []
            else:
                for dp in OEGetSDDataIter(mol1):
                    OESetSDData(molDict[smiles],dp.GetTag(),dp.GetValue())
            if vendor_name not in vendorDict[smiles]:
                vendorDict[smiles].append(vendor_name)
        ifs.close()

    ofs = oemolostream()
    ofs.open("output_all.sdf")
    for smiles in molDict:
        mol = molDict[smiles]
        OESetSDData(mol,"Vendors","|".join(vendorDict[smiles]))
        if bioNumberDict.has_key(smiles):
            OESetSDData(mol,"BIO-NUMBER",bioNumberDict[smiles])
        OEWriteMolecule(ofs,mol)
    ofs.close()



