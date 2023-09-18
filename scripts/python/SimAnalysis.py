#!/usr/bin/env python
import sys
from openeye.oechem import *
from openeye.oegraphsim import *

if len(sys.argv) != 4:
    OEThrow.Usage("%s <database> <input> <output>" % sys.argv[0])


# load molecules

probes = []
mol = OEGraphMol()
ifs = oemolistream()
ifs.open(sys.argv[2])
while OEReadMolecule(ifs,mol):
    probes.append(OEGraphMol(mol))
ifs.close()

ifs = oemolistream()
if not ifs.open(sys.argv[1]):
    OEThrow.Fatal("Cannot open database molecule file!")
moldb = OEMolDatabase(ifs)
nrmols = moldb.GetMaxMolIdx()

# generate fingerprint

fpdb = OEFPDatabase(OEFPType_Path)

emptyfp = OEFingerPrint()
emptyfp.SetFPTypeBase(fpdb.GetFPTypeBase())

for idx in range(0, nrmols):
    mol = OEGraphMol()
    if moldb.GetMolecule(mol, idx):
        fpdb.AddFP(mol)
    else:
        fpdb.AddFP(emptyfp)

nrfps = fpdb.NumFingerPrints()

#reading probes
ofs = oemolostream()
ofs.SetFormat(OEFormat_SDF)
ofs.open(sys.argv[3])
print len(probes)

for idx,probe in enumerate(probes):
    hits = fpdb.GetSortedScores(probe, 1)
    for hit in hits:
        hitMol = OEGraphMol()
        result = moldb.GetMolecule(hitMol, hit.GetIdx())
        hit_smiles = OEMolToSmiles(hitMol)
        if result:
            print idx
            OESetSDData(probe,"reference_name",hitMol.GetTitle())
            OESetSDData(probe,"reference_smiles",hit_smiles)
            OESetSDData(probe,"tanimoto_score","%5.2f"%hit.GetScore())
            OEWriteMolecule(ofs,probe)
        else:
            print "failed"
            print hit_smiles
ofs.close()
