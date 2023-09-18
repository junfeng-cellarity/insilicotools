#!/usr/bin/env python

import rdkit
from rdkit import Chem
# from FPSim2.io import create_db_file
# create_db_file('/home/jfeng/Database/OnSiteScreenDb/skyhawk_chemdiv.sdf', 'skyhawk_chemdiv.h5', 'Morgan', {'radius': 3, 'nBits': 2048}, gen_ids=True)
import sys
from FPSim2 import FPSim2Engine
database = "/home/jfeng/Database/skyhawk_chemdiv.h5"

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print("usage:%s input.sdf output.sdf cutoff"%sys.argv[0])
        print("remove compounds in input.sdf that has nearest neighbor (tanimoto > cutoff) in RGX database, save result to output.sdf")
    else:
        cutoff = float(sys.argv[3])
        fpe = FPSim2Engine(database, in_memory_fps=False)
        sd_reader = Chem.SDMolSupplier(sys.argv[1])
        sd_writer = Chem.SDWriter(sys.argv[2])
        for mol in sd_reader:
            try:
                query = Chem.MolToSmiles(mol)
                results = fpe.on_disk_similarity(query, cutoff, n_workers=25)
                if len(results) == 0:
                    sd_writer.write(mol)
            except:
                continue
        sd_writer.close()
