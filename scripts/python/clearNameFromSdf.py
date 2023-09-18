__author__ = 'jfeng1'

from openeye.oechem import *

ifs = oemolistream()
ifs.open("chembl_test.sdf")
ofs = oemolostream()
ofs.SetFormat(OEFormat_SDF)
ofs.open("output.sdf")

mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    mol.SetTitle("")
    OEWriteMolecule(ofs,mol)
ifs.close()
ofs.close()