import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
#drugs = [ mol for mol in Chem.SDMolSupplier( "/home/jfeng/Database/AldrichMarketSelect/analysis2.sdf" ) if mol != None ]
# colorDict = {
#     "RGX":"red",
#     "Library":"blue",
#     "Diverse":"green"
# }
drugs = [ mol for mol in Chem.SDMolSupplier( "/home/jfeng/Downloads/tsne_0906_2.sdf" ) if mol != None ]
colorDict = {
    "EnamineLibrary":"green",
    "Skyhawk":"green",
    "AllEnumerated":"green",
    "Library_2019_2":"yellow",
    "MarketSelectDiverse50k":"blue",
    "Risdiplam":"red",
    "Branaplam":"red"
}

orderDict = {
    "EnamineLibrary":2,
    "Skyhawk":3,
    "AllEnumerated":1,
    "MarketSelectDiverse50k":0,
    "Risdiplam":4,
    "Branaplam":5,
    "Library_2019_2":2
}


def calc_fp_arr( mols ):
    source = []
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 3 )
        color_key = mol.GetProp("Source")
        #color_key = mol.GetProp("class")
        source.append(colorDict[color_key])
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist ),source

res,source =calc_fp_arr(drugs)
model = TSNE( n_components=2, init="pca",  n_iter=2000,  )
#model = TSNE( n_components=2,random_state=42 ,  n_iter=1000, metric="jaccard")
x = model.fit_transform( res )
x = pd.DataFrame( x )
x.columns = ['A1', 'A2']

sd_writer = Chem.SDWriter("/home/jfeng/tsne_tmp_2.sdf")
for id,mol in enumerate(drugs):
    order_key = mol.GetProp("Source")
    order = orderDict[order_key]
    mol.SetProp("Order",str(order))
    mol.SetProp("X",str(x['A1'][id]))
    mol.SetProp("Y",str(x['A2'][id]))
    sd_writer.write(mol)
sd_writer.close()

plt.scatter(x['A1'],x['A2'],color=source)
plt.show()