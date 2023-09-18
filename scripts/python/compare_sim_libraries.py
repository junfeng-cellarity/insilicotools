#!/usr/bin/env python

import rdkit
import sys
import rdkit.Chem
from matplotlib import pyplot
import numpy as np
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: Compare library similarity histogram")
        print("Usage: %s input1.sdf input2.sdf")
    else:
        sdf_reader_1 = rdkit.Chem.SDMolSupplier(sys.argv[1])
        sim_1 = []
        for mol in sdf_reader_1:
            tanimoto = float(mol.GetProp("Tanimoto"))
            sim_1.append(tanimoto)

        sdf_reader_2 = rdkit.Chem.SDMolSupplier(sys.argv[2])
        sim_2 = []
        for mol in sdf_reader_2:
            tanimoto = float(mol.GetProp("Tanimoto"))
            sim_2.append(tanimoto)

        bins = np.linspace(0,1.0,20)
        pyplot.hist(sim_1,bins,alpha=0.5,label="targetmol",color="red")
        pyplot.hist(sim_2,bins,alpha=0.5,label="custom",color="green")
        pyplot.legend(loc="upper right")
        pyplot.show()
