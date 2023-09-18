#!/usr/bin/env python
import os
import rdkit
from rdkit import Chem
from FPSim2.io import create_db_file
import sys
from FPSim2 import FPSim2Engine

if __name__ == "__main__":
    if len(sys.argv)!=5:
        print("usage:%s input.sdf datbase.sdf output.sdf cutoff"%sys.argv[0])
        print("remove compounds in input.sdf that has nearest neighbor (tanimoto > cutoff) in database.sdf, save result to output.sdf")
    else:
        database_sdf = sys.argv[2]
        fingerprint_file = "%s.h5"%database_sdf.rsplit(".",1)[0]
        create_db_file(database_sdf, fingerprint_file, 'Morgan', {'radius': 3, 'nBits': 2048}, gen_ids=True)
        cutoff = float(sys.argv[4])
        fpe = FPSim2Engine(fingerprint_file, in_memory_fps=False)
        sd_reader = Chem.SDMolSupplier(sys.argv[1])
        sd_writer = Chem.SDWriter(sys.argv[3])
        for mol in sd_reader:
            try:
                query = Chem.MolToSmiles(mol)
                results = fpe.on_disk_similarity(query, cutoff, n_workers=25)
                if len(results) == 0:
                    sd_writer.write(mol)
            except:
                continue
        sd_writer.close()
