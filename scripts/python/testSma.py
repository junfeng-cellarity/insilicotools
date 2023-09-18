#!/usr/bin/env python
from openeye.oechem import *
#smiles = "n1c2c(nc[nH]2)c([nH]c1N)=S"
sdf = """
  Mrv15c1403021615122D

 11 12  0  0  0  0            999 V2000
    0.7145   -0.4125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0000    0.8250    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7846    1.0799    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2695    0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7846   -0.2549    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    1.2375    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    0.8250    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4289    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1434   -0.4125    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.7145    2.0625    0.0000 S   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  4  0  0  0  0
  2  3  4  0  0  0  0
  3  4  4  0  0  0  0
  4  5  4  0  0  0  0
  5  6  4  0  0  0  0
  2  6  4  0  0  0  0
  3  7  4  0  0  0  0
  7  8  4  0  0  0  0
  8  9  4  0  0  0  0
  1  9  4  0  0  0  0
  9 10  1  0  0  0  0
  7 11  2  0  0  0  0
M  STY  2   1 DAT   2 DAT
M  SAL   1  1   6
M  SDT   1 MRV_IMPLICIT_H
M  SDD   1     0.0000    0.0000    DR    ALL  0       0
M  SED   1 IMPL_H1
M  SAL   2  1   8
M  SDT   2 MRV_IMPLICIT_H
M  SDD   2     0.0000    0.0000    DR    ALL  0       0
M  SED   2 IMPL_H1
M  END
"""

smiles = "NC1=NC2=C(N=CN2)C(=S)N1"
smiles = "C1=NC2=NC(=NC(=S)C2=N1)N"
smiles = "Nc1nc2[nH]cnc2c(=S)[nH]1"
smarts = "[N](C)=C"
subsearch = OESubSearch(smarts)

mol = OEGraphMol()
ifs = oemolistream()
ifs.SetFormat(OEFormat_SDF)
ifs.openstring(sdf)
OEReadMolecule(ifs,mol)
#print("Number of aromatic atoms =", OECount(mol, OEIsAromaticAtom()))
OEAssignAromaticFlags(mol)
print OEMolToSmiles(mol)
print("Number of aromatic atoms =", OECount(mol, OEIsAromaticAtom()))

OEPrepareSearch(mol,subsearch)

m = subsearch.Match(mol,True)
for idx,match in enumerate(m):
    print idx
