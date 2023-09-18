__author__ = 'jfeng1'
from openeye.oechem import *
import glob
import os

ofs = oemolostream()
ofs.openstring()
ofs.SetFormat(OEFormat_SDF)

sdfs = glob.glob("*.sdf")
for sdf in sdfs:
    reagent_class = os.path.splitext(os.path.basename(sdf))[0]
    ifs = oemolistream()
    ifs.open(sdf)
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        OESetSDData(mol,"Reagent Class",reagent_class)
        OEWriteMolecule(ofs,mol)
print ofs.GetString()
