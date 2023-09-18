#!/usr/bin/env python
from openeye.oechem import *
import sys,os
from sqlitedict import SqliteDict

if __name__=="__main__":
    if len(sys.argv)!=3:
        print "Add Evotec Inventory Info"
        print "Usage:%s input.sdf output.sdf"%sys.argv[0]
    else:

        inputfile = sys.argv[1]
        outputfile = sys.argv[2]

        if os.path.isfile(inputfile):
            db = SqliteDict("/Users/jfeng1/BiogenDB/Evotec/EvotecInventory.db",autocommit=True)
            ifs = oemolistream()
            ifs.open(inputfile)
            mol = OEGraphMol()
            ofs = oemolostream()
            ofs.open(outputfile)
            while OEReadMolecule(ifs,mol):
                #bio_number = OEGetSDData(mol,"BIO-NUMBER")
                bio_number = mol.GetTitle()
                if db.has_key(bio_number):
                    OESetSDData(mol,"has_inventory","Yes")
                    OESetSDData(mol,"evotec_inventory",db[bio_number])
                else:
                    OESetSDData(mol,"has_inventory","No")
                    OESetSDData(mol,"evotec_inventory","None")
                OEWriteMolecule(ofs,mol)
            ofs.close()
            ifs.close()
