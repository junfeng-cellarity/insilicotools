#!/usr/bin/env python
from openeye.oechem import *
import numpy as np
import sys

def calc_moments_of_inertia(mol):
    center_of_mass = OEFloatArray(3)
    OEGetCenterOfMass(mol, center_of_mass, True)
    print center_of_mass
    inertia = np.zeros((3, 3))

    for atm in mol.GetAtoms():
        cmx, cmy, cmz = np.array(mol.GetCoords(atm)) - np.array(center_of_mass)
        elem = atm.GetAtomicNum()
        mass = OEGetAverageWeight(elem)
        inertia += - mass * np.matrix(
            [[0, -cmz, cmy],
             [cmz, 0, -cmx],
             [-cmy, cmx, 0]]) ** 2
        print inertia

    return np.linalg.eigvalsh(inertia)

def moi_coords(mol):
    inertia = calc_moments_of_inertia(mol)
    print inertia
    moi = sorted(inertia)
    if moi[2] > 0.0:
        res = moi[0:2]/moi[2]
    else:
        res = [0.0,0.0]
    return res

if __name__=="__main__":
    if len(sys.argv)!=2:
        print "Usage:%s input.sdf"%(sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            print moi_coords(mol)
