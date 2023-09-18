#!/usr/bin/env python
import os,sys
import csv
from openeye.oechem import *
if __name__ == "__main__":
    csvfile = open("/home/jfeng/Downloads/salts.csv","r")
    saltDict = {}
    csvreader = csv.DictReader(csvfile)
    for row in csvreader:
        mol = OEGraphMol()
        smiles = row['SMILES']
        OEParseSmiles(mol,smiles)
        mw = OECalculateMolecularWeight(mol)
        salt = row['Name']
        saltDict[salt]=mw

    sdfiles = ['OTAVA_CDD_Output.sdf','Life Chemical Library_CDD_output.sdf']
    outputfiles = ["OTAVA_CDD_new.sdf","Life_Chemical_New.sdf"]
    for id,f in enumerate(sdfiles):
        f = os.path.join("/home/jfeng/Downloads",f)
        of = os.path.join("/home/jfeng/",outputfiles[id])
        ofs = oemolostream()
        ofs.open(of)
        ifs = oemolistream()
        ifs.open(f)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            if id==0:
                barcode = str(int(float(OEGetSDData(mol,"Barcode"))))
                OESetSDData(mol,"Barcode",barcode)
            parent_mw = float(OEGetSDData(mol,"Molecular weight (g/mol)"))
            salt_name = OEGetSDData(mol,"Batch Salt")
            if not salt_name.startswith("No Salt"):
                #print(salt_name)
                eq,name = salt_name.split(" ")
                eq = float(eq)
                name = name.strip()
                #print(name)
                if name in saltDict:
                    mw = parent_mw + eq*saltDict[name]
                    OESetSDData(mol,"Total MW","%f"%mw)
            else:
                OESetSDData(mol,"Total MW","%f"%parent_mw)
            OEWriteMolecule(ofs,mol)
        ofs.close()

