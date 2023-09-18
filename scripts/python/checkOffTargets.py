__author__ = 'jfeng1'
import sys
from openeye.oechem import *

ifs = oemolistream()
ifs.open("/Users/jfeng1/two-sets-100nMa.sdf")
mol = OEGraphMol()
n = 0
molNames = []

while OEReadMolecule(ifs,mol):
    numOfPrimary = 0
    numOfOffTargets = 0
    for dp in OEGetSDDataIter(mol):
        tag  = dp.GetTag()
        if tag == "Test-conc":
           continue
        try:
            value = int(dp.GetValue())
            if value < 5:
                numOfPrimary += 1
            elif value < 20:
                numOfOffTargets += 1
        except:
            continue

    if numOfPrimary >= 1:
        if numOfPrimary+numOfOffTargets-1 <= 5:
            if mol.GetTitle() not in molNames:
                molNames.append(mol.GetTitle())

print len(molNames)
print molNames




