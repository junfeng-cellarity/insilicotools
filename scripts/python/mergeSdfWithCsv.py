#!/usr/bin/env python
from openeye.oechem import *
import sys
import csv
if __name__=="__main__":
    if len(sys.argv)!=4:
        print("Usage:%s input.sdf input.csv output.sdf"%sys.argv[0])
        print("merge sdf with input.csv, first line of csv file contains tag names")
    else:
        csvfile = open(sys.argv[2],"r")
        reader = csv.reader(csvfile)
        id_name = None
        tag_name = None
        dict = {}
        for row_id,row in enumerate(reader):
            if row_id ==0:
                id_name = row[0]
                tag_name = row[1]
            else:
                dict[row[0]] = row[1]
        csvfile.close()

        ifs = oemolistream()
        ifs.open(sys.argv[1])
        ofs = oemolostream()
        ofs.open(sys.argv[3])
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            if mol.GetTitle() not in dict:
                continue
            OESetSDData(mol,tag_name,dict[mol.GetTitle()])
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()