#!/usr/bin/env python
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import sys

if len(sys.argv)!=3:
    print("Usage:%s input.sdf output.sdf"%sys.argv[0])
    sys.exit(1)

drugs = [ mol for mol in Chem.SDMolSupplier( sys.argv[1] ) if mol != None ]

def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 3 )
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist )

res =calc_fp_arr(drugs)
model = TSNE( n_components=2, init="pca",  n_iter=2000 )
#model = TSNE( n_components=2,random_state=42 ,  n_iter=1000 )
x = model.fit_transform( res )
x = pd.DataFrame( x )
x.columns = ['A1', 'A2']

sd_writer = Chem.SDWriter(sys.argv[2])
for id,mol in enumerate(drugs):
    mol.SetProp("X",str(x['A1'][id]))
    mol.SetProp("Y",str(x['A2'][id]))
    sd_writer.write(mol)
sd_writer.close()

# plt.scatter(x['A1'],x['A2'])
# plt.show()