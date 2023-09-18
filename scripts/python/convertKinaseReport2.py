#!/usr/bin/env python
__author__ = 'jfeng1'

from openeye.oechem import *
ignored_tags = ["Name","Concat Distinct;Project","BIO Number"]
sdf1000 = "/Users/jfeng1/Datasets/forJaclyn/BID083-01-p-00003_Percent_Control_1000nM_96.sdf"
sdf100 = "/Users/jfeng1/Datasets/forJaclyn/BID083-01-p-00003_Percent_Control_100nM_96.sdf"

ifs = oemolistream()
ifs.open(sdf1000)
mol = OEGraphMol()
molList = []
kinaseTags = []
molDict = {}
dataDict = {}
while OEReadMolecule(ifs,mol):
    bio_number = mol.GetTitle()
    molDict[bio_number] = OEGraphMol(mol)
    molList.append(bio_number)
    for dp in OEGetSDDataIter(mol):
        if not dp.GetTag() in ignored_tags:
            value = dp.GetValue()
            data_tag = "%s_%s_%d"%(bio_number,dp.GetTag(),1000)
            dataDict[data_tag] = value
            if not dp.GetTag() in kinaseTags:
                kinaseTags.append(dp.GetTag())
ifs.close()

ifs.open(sdf100)
while OEReadMolecule(ifs,mol):
    bio_number = mol.GetTitle()
    for dp in OEGetSDDataIter(mol):
        if not dp.GetTag() in ignored_tags:
            value = dp.GetValue()
            data_tag = "%s_%s_%d"%(bio_number,dp.GetTag(),100)
            dataDict[data_tag] = value
ifs.close()

ofs = oemolostream()
ofs.SetFormat(OEFormat_SDF)
ofs.open("kinase4.sdf")
for bio_number in molList:
    if not molDict.has_key(bio_number):
        continue

    mol = molDict[bio_number]
    numHits = 0
    hasProlem = 0
    hasActivity = 0
    for kinase_name in kinaseTags:
        OEDeleteSDData(mol,kinase_name)
        concentration1 = 100
        data_key_1 = "%s_%s_%d"%(bio_number,kinase_name,concentration1)
        value1 = dataDict[data_key_1]

        concentration2 = 1000
        data_key_2 = "%s_%s_%d"%(bio_number,kinase_name,concentration2)
        value2 = dataDict[data_key_2]

        v1 = float(value1)
        v2 = float(value2)
        if v1 < 30 or v2 < 30:
            if v1 < 30 and v2 > 30:
                hasProlem = 1
            else:
                OESetSDData(mol,"%s_%s"%(kinase_name,concentration1),value1)
                OESetSDData(mol,"%s_%s"%(kinase_name,concentration2),value2)
                if v1<30 and v2<30:
                    dose_response = 1
                numHits += 1
                hasActivity = 1
    if hasActivity:
        OESetSDData(mol,"Number of Hits","%d"%numHits)
        OESetSDData(mol,"Has Problem","%d"%hasProlem)
        OEWriteMolecule(ofs,mol)
ofs.close()
