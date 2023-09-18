import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from MulticoreTSNE import MulticoreTSNE as TSNE
import matplotlib.pyplot as plt
drugs = [ mol for mol in Chem.SDMolSupplier( "/home/jfeng/Database/AldrichMarketSelect/analysis2.sdf" ) if mol != None ]
colorDict = {
    "RGX":"red",
    "Library":"blue",
    "Diverse":"green"
}
zorderDict = {
    "RGX":10,
    "Library":5,
    "Diverse":1
}

zorder = []
def calc_fp_arr( mols ):
    source = []
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 3 )
        color_key = mol.GetProp("Source")
        source.append(colorDict[color_key])
        zorder.append(zorderDict[color_key])
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist ),source

res,source =calc_fp_arr(drugs)
model = TSNE(n_jobs=80,n_components=2,n_iter=2000)
#model = TSNE( n_components=2, init="pca",  n_iter=2000 )
#model = TSNE( n_components=2,random_state=42 ,  n_iter=1000 )
x = model.fit_transform( res )
x = pd.DataFrame( x )
x.columns = ['A1', 'A2']

sd_writer = Chem.SDWriter("/home/jfeng/tsne.sdf")
for id,mol in enumerate(drugs):
    mol.SetProp("X","%f"%x['A1'][id])
    mol.SetProp("Y","%f"%x['A2'][id])
    sd_writer.write(mol)
sd_writer.close()

plt.scatter(x['A1'],x['A2'],color=source,alpha=0.1)
plt.show()