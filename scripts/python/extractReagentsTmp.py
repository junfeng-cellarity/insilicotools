#!/usr/bin/env python

from openeye.oechem import *
import sys,os
import glob
from sqlitedict import SqliteDict
reagentDirectory = "/Users/jfeng1/WuXi_Library/"

if __name__=="__main__":
    # reagent_db = SqliteDict("/Users/jfeng1/Databases/ReagentUsage.db",autocommit=True)
    # reagent_db.clear()
    reagent_db = {}
    files = glob.glob("%s/*.sdf"%reagentDirectory)
    for f in files:
        ifs = oemolistream()
        ifs.open(f)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            r1,r2 = mol.GetTitle().split("_")
            if reagent_db.has_key(r1):
                reagent_db[r1] = reagent_db[r1]+1
            else:
                reagent_db[r1] = 1

            if reagent_db.has_key(r2):
                reagent_db[r2] = reagent_db[r2]+1
            else:
                reagent_db[r2] = 1
        ifs.close()
    reagents = []
    for key in reagent_db.keys():
        reagents.append([key, reagent_db[key]])
    #sorted(plates,key=lambda plate:plate[1])
    reagents = sorted(reagents,key=lambda reagent:reagent[1])
    for r in reagents:
        print r




