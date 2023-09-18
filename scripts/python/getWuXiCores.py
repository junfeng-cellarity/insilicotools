#!/usr/bin/env python
from openeye.oechem import *

mol = OEGraphMol()

coreSmiDict20g = {}
ifs = oemolistream()
ifs.open("/Users/jfeng1/Datasets/WuXi/WuXi Template list (more than 20g)_382cpds.sdf")
while OEReadMolecule(ifs,mol):
    coreId = OEGetSDData(mol,"WUXIID")
    coreSmiDict20g[coreId] = OEMolToSmiles(mol)
ifs.close()

coreSmiDict5g = {}
ifs = oemolistream()
ifs.open("/Users/jfeng1/Datasets/WuXi/WuXiCore_5g_Nov26_2014.sdf")
while OEReadMolecule(ifs,mol):
    coreId = OEGetSDData(mol,"Text3")
    coreSmiDict5g[coreId] = OEMolToSmiles(mol)
ifs.close()

ifs = oemolistream()
ifs.open("/Users/jfeng1/Datasets/WuXi/EnumeratedLibrary/Set1/WuXi_2600_Set_1-3.sdf")


ofs = oemolostream()
ofs.open("/Users/jfeng1/Datasets/WuXi/EnumeratedLibrary/Set1/WuXi_2600_Set_1-3_With_Core.sdf")
coreList = []
while OEReadMolecule(ifs,mol):
    coreName = OEGetSDData(mol,"Alias_ID").strip().split("_")[0]
    if coreSmiDict5g.has_key(coreName):
        OESetSDData(mol,"coreSmiles",coreSmiDict5g[coreName])
    elif coreSmiDict20g.has_key(coreName):
        OESetSDData(mol,"coreSmiles",coreSmiDict20g[coreName])
    if coreName not in coreList:
        coreList.append(coreName)
    OEWriteMolecule(ofs,mol)
ofs.close()

for core in coreList:
    if coreSmiDict5g.has_key(core):
        print coreSmiDict5g[core], core, "5g", 1
        continue
    if coreSmiDict20g.has_key(core):
        print coreSmiDict20g[core], core, "20g", 1
        continue
    else:
        print "Not Found for %s"%core


ifs.close()
