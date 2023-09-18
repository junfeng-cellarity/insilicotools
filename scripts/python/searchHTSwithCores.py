#!/usr/bin/env python
from openeye.oechem import *
import os,sys

if len(sys.argv)<3:
    print "Usage:%s coreToSearch.sdf output.sdf"%sys.argv[0]
    sys.exit(1)

# HTS_SDF = "/Users/jfeng1/Databases/HTS/hts.sdf"
HTS_SDF = "/Users/jfeng1/Datasets/CoresCollection/cores.sdf"
CORE_SDF = sys.argv[1]

ifs = oemolistream()
ifs.open(CORE_SDF)
coreMol = OEGraphMol()
coreNames = []
coreSmilesDict = {}
coreDict = {}
result = {}

while OEReadMolecule(ifs,coreMol):
    coreName = OEGetSDData(coreMol,"Scaffold ID")
    coreNames.append(coreName)
    result[coreName] = []
    subsearch = OESubSearch(coreMol,OEExprOpts_DefaultAtoms,OEExprOpts_DefaultBonds)
    coreDict[coreName] = subsearch
    coreSmilesDict[coreName] = OEMolToSmiles(coreMol)
ifs.close()

ofs = oemolostream()
ofs.open(sys.argv[2])

ifs = oemolistream()
ifs.open(HTS_SDF)
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    for coreName in coreNames:
        subsearch = coreDict[coreName]
        OEPrepareSearch(mol,subsearch)
        OEAddExplicitHydrogens(mol)
        if subsearch.SingleMatch(mol):
            OESetSDData(mol,"MatchedCore",coreName)
            OESetSDData(mol,"MatchedCoreSmiles",coreSmilesDict[coreName])
            OEWriteMolecule(ofs,mol)
            result[coreName].append(mol.GetTitle())
ifs.close()
ofs.close()

for key in result:
    if len(result[key])>0:
        print key,len(result[key])

