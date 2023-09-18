#!/usr/bin/env python

from openeye.oechem import *
from openeye.oemolprop import *
import sys, os, xmlrpclib, glob, random,operator
from openeye.oegraphsim import *
from sqlitedict import SqliteDict
from progressbar import ProgressBar
import random


def extractPlatesFromFile(plateFile):
    lines = open(plateFile,"r").read().splitlines()
    plates = []
    for line in lines:
        args = line.split()
        plates.append(args[0])
    return plates

if __name__=="__main__":
    ifs = oemolistream()
    ifs.open("/Users/jfeng1/Datasets/OGA/softfocus_diverse_uniq.sdf.gz")
    smiDict = {}
    plateDict = {}
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        smiles = OEMolToSmiles(mol)
        name = OEGetSDData(mol,"Abase_ID")
        plate = OEGetSDData(mol,"plate_ID")
        smiDict[smiles] = name
        plateDict[name] = plate
    ifs.close()

    ifs = oemolistream()
    ifs.open("/Users/jfeng1/Datasets/Glycomimetic/glycomimetic library.sdf")
    glyDict = {}
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        smiles = OEMolToSmiles(mol)
        name = OEGetSDData(mol,"IDNUMBER")
        glyDict[smiles] = name
    ifs.close()

    biogenDict = {}
    ifs = oemolistream()
    ifs.open("/Users/jfeng1/Datasets/OGA/evotec_biib_set.sdf")
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        smiles = OEMolToSmiles(mol)
        name = OEGetSDData(mol,"BIO Number")
        biogenDict[smiles] = name
    ifs.close()

    rocs_plates = extractPlatesFromFile("/Users/jfeng1/Datasets/OGA/rocs_plates.txt")
    leadlike_plates = extractPlatesFromFile("/Users/jfeng1/Datasets/OGA/lead_like_plates.txt")

    ofs = oemolostream()
    ofs.open("/Users/jfeng1/Datasets/OGA/OGA_dose_responses.sdf")

    ifs = oemolistream()
    ifs.open("/Users/jfeng1/Datasets/OGA/OGA_dose_response_1728cmpds.sdf")
    mol = OEGraphMol()
    numRocs = 0
    numLeadLike = 0
    while OEReadMolecule(ifs,mol):
        name = OEGetSDData(mol,"Corporate ID")
        smiles = OEMolToSmiles(mol)
        if plateDict.has_key(name):
            plateId = plateDict[name]
            if plateId in rocs_plates:
                OESetSDData(mol,"Plate_Source","ROCS_PLATES")
                numRocs += 1
            elif plateId in leadlike_plates:
                OESetSDData(mol,"Plate_Source","LEADLIKE_PLATES")
                numLeadLike += 1
            else:
                OESetSDData(mol,"Plate_Source","CRL_OTHER")
        elif glyDict.has_key(smiles):
            OESetSDData(mol,"Plate_Source","BIOGEN_HIS(Glycomimetic)")

        elif smiDict.has_key(smiles):
            name = smiDict[smiles]
            if plateDict.has_key(name):
                plateId = plateDict[name]
                if plateId in rocs_plates:
                    OESetSDData(mol,"Plate_Source","ROCS_PLATES")
                    numRocs += 1
                elif plateId in leadlike_plates:
                    OESetSDData(mol,"Plate_Source","LEADLIKE_PLATES")
                    numLeadLike += 1
                else:
                    OESetSDData(mol,"Plate_Source","CRL_OTHER")
            else:
                OESetSDData(mol,"Plate_Source","Unknown")
        else:
            if biogenDict.has_key(smiles):
                OESetSDData(mol,"Plate_Source","BIOGEN_HTS")
            else:
                OESetSDData(mol,"Plate_Source","Unknown")

        OEWriteMolecule(ofs,mol)
    ofs.close()

    print numRocs,numLeadLike
