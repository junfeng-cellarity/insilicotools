#!/usr/bin/env python
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
import progressbar

if len(sys.argv) != 4:
    print("%s <database> <input> <output>" % sys.argv[0])
    sys.exit(1)

# load molecules
probes = [mol for mol in Chem.SDMolSupplier(sys.argv[2]) if mol != None]
database = [mol for mol in Chem.SDMolSupplier(sys.argv[1]) if mol != None]

# generate fingerprint
probe_dict = {}
fp_probes = []
for idx,mol in enumerate(probes):
    fp = Chem.AllChem.GetMorganFingerprint(mol, 3)
    fp_probes.append(fp)
    probe_dict[idx] = mol

db_dict = {}
fp_db = []
for idx,mol in enumerate(database):
    fp = Chem.AllChem.GetMorganFingerprint(mol,3)
    fp_db.append(fp)
    db_dict[idx] = mol

#reading probes
writer = Chem.SDWriter(sys.argv[3])
progressbar = progressbar.ProgressBar(len(fp_probes))
progressbar.start()
for probe_idx,probe_fp in enumerate(fp_probes):
    progressbar.update(probe_idx)
    best_tanimoto = 0
    best_idx = -1
    for db_idx,db_fp in enumerate(fp_db):
        tanimoto = DataStructs.TanimotoSimilarity(probe_fp,db_fp)
        if tanimoto>best_tanimoto:
            best_tanimoto = tanimoto
            best_idx = db_idx
    mol = probe_dict[probe_idx]
    mol.SetProp("Tanimoto","%5.2f"%best_tanimoto)
    if best_idx in db_dict:
        mol.SetProp("SimMolName",db_dict[best_idx].GetProp("_Name"))
        mol.SetProp("SimSmiles",Chem.MolToSmiles(db_dict[best_idx]))
    else:
        print("Not found %d"%best_idx)
    writer.write(mol)
