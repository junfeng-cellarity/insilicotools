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
    rocs_plates = extractPlatesFromFile("/Users/jfeng1/Datasets/OGA/rocs_plates.txt")
    leadlike_plates = extractPlatesFromFile("/Users/jfeng1/Datasets/OGA/lead_like_plates.txt")

    ifs = oemolistream()
    ifs.open("/Users/jfeng1/Datasets/OGA/softfocus_diverse_uniq.sdf.gz")
    nRocs = 0
    nLeadLike = 0
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        plate = OEGetSDData(mol,"plate_ID")
        if plate in rocs_plates:
            nRocs += 1
        if plate in leadlike_plates:
            nLeadLike +=1
    print nRocs,nLeadLike
    ifs.close()

