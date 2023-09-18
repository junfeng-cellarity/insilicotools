#!/usr/bin/env python
from openeye.oechem import *
import glob,os

molFiles = glob.glob("*.mol")

ofs = oemolostream()
ofs.openstring()
ofs.SetFormat(OEFormat_SDF)
for molfile in molFiles:
    molName = os.path.splitext(molfile)[0]
    mol = OEGraphMol()
    ifs = oemolistream()
    ifs.SetFormat(OEFormat_MDL)
    ifs.open(molfile)
    OEReadMolecule(ifs,mol)
    mol.SetTitle(molName)
    OEWriteMolecule(ofs,mol)
    ifs.close()
print ofs.GetString()
ofs.close()


