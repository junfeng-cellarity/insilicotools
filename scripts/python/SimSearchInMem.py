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
existing_hits = []
ofs = oemolostream()
ofs.SetFormat(OEFormat_SDF)
ofs.open(sys.argv[3])
for probe in probes:
    hits = fpdb.GetSortedScores(probe, 10)
    for hit in hits:
        hitMol = OEGraphMol()
        result = moldb.GetMolecule(hitMol, hit.GetIdx())
        if result:
            hit_smiles = OEMolToSmiles(hitMol)
            if hit_smiles not in existing_hits:
                existing_hits.append(hit_smiles)
                OESetSDData(hitMol,"reference_name",probe.GetTitle())
                OESetSDData(hitMol,"reference_smiles",OEMolToSmiles(probe))
                OESetSDData(hitMol,"tanimoto_score","%5.2f"%hit.GetScore())
                OEWriteMolecule(ofs,hitMol)
        else:
            print "failed"
ofs.close()
