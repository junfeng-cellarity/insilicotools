#!/usr/bin/env python

from openeye.oechem import *
import sys,os

if len(sys.argv)!=3:
    print ("%s input.sdf output.sdf"%sys.argv[0])
    sys.exit(1)

ifs = oemolistream()
lowestEnergyDict = {}
ifs.open(sys.argv[1])
mol = OEGraphMol()
while(OEReadMolecule(ifs,mol)):
    title = mol.GetTitle()
    energy = float(OEGetSDData(mol,"mmff94smod_NoEstat"))
    if title not in lowestEnergyDict:
        lowestEnergyDict[title] = energy
    else:
        if energy < lowestEnergyDict[title]:
            lowestEnergyDict[title] = energy
ifs.close()

ifs.open(sys.argv[1])
ofs = oemolostream()
ofs.open(sys.argv[2])
mol = OEGraphMol()
while(OEReadMolecule(ifs,mol)):
    title = mol.GetTitle()
    energy = float(OEGetSDData(mol,"mmff94smod_NoEstat"))
    lowest_energy = lowestEnergyDict[title]
    OESetSDData(mol,"Relative_Energy","%5.2f"%(energy-lowest_energy))
    OEWriteMolecule(ofs,mol)
ifs.close()
ofs.close()