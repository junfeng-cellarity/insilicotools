#!/usr/bin/env python

controlFile = open("/Users/jfeng1/Datasets/Eric/HbF_all_controls.txt","r")
sodium_butyrate_smiles = "CCCC(=O)[O-].[Na+]"
sodium_butyrate_formula = "C4H7NaO2"
sodium_butyrate_mw = 110.09
dmso_smiles = "CS(=O)C"
dmso_formula = "C2H6OS"
dmso_mw = 78.13

lines = controlFile.readlines()
for lineno,line in enumerate(lines):
    if lineno ==0:
        continue
    else:
        plateId,well,positiveControl,negativeControl = line.strip().split("\t")
        print "%s\t%s\t\t\t%f\t%s\t%s\t%s1"%(sodium_butyrate_smiles,positiveControl,sodium_butyrate_mw,sodium_butyrate_formula,plateId,well)
        print "%s\t%s\t\t\t%f\t%s\t%s\t%s2"%(dmso_smiles,negativeControl,dmso_mw,dmso_formula,plateId,well)
controlFile.close()
