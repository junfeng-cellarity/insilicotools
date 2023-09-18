#!/usr/bin/env python
import pandas, csv
file = open("/home/jfeng/Downloads/scaffold_all_genes_expression_expected_count.csv", "r")
dict_reader = csv.DictReader(file, dialect="excel-tab")
gene_list = []
compound_list = []
batch_list = []
experiment_list = []
position_list = []
gene_expressions = []
id = 0
uniq_gene_list = []
for row in dict_reader:
    uniq_gene_list.append((row['gene']))
    columns = row.keys()
    for column in columns:
        if column != "gene":
            batch,position,experiment,rgx_number,exp_id = column.split("_")
            gene_list.append(row['gene'])
            compound_list.append(rgx_number)
            batch_list.append(batch)
            experiment_list.append(experiment)
            gene_expressions.append(float(row[column]))
            position_list.append(position)
print(len(gene_list),len(gene_expressions),len(position_list),len(batch_list),len(compound_list))

df = pandas.DataFrame.from_dict({"batch":batch_list,"gene":gene_list,"molecular_name":compound_list,"gene_expressions":gene_expressions,"exp_ids":experiment_list})

#df_1 = df.query("batch=='Batch2'")
df_1 = df
df_1 = df_1.filter(items=["molecular_name","gene","gene_expressions","batch"]).groupby(["molecular_name","batch","gene"]).mean()
print(df_1)
control = {}
data = []
uniq_compound_list = []

for index,row in df_1.iterrows():
    compound_name = index[0]
    batch_name = index[1]
    gene_name = index[2]
    if compound_name not in uniq_compound_list:
        uniq_compound_list.append(compound_name)
    if index[0]=="DMSO":
        if batch_name not in control:
            control[batch_name] = {}
        control[batch_name][gene_name] = row[0]

for index,row in df_1.iterrows():
    if index[0] not in uniq_compound_list:
        uniq_compound_list.append(index[0])

data = []
for c in uniq_compound_list:
    data.append({})

# import math
for index,row in df_1.iterrows():
    compound_name = index[0]
    batch_name = index[1]
    gene_name = index[2]
    idx = uniq_compound_list.index(compound_name)
    dmso = control[batch_name][gene_name]
    value = row[0]-dmso
    data[idx][gene_name] = value
new_data = []
for idx in range(len(data)):
    cmpd_data = []
    for gene in uniq_gene_list:
        cmpd_data.append(data[idx][gene])
    new_data.append(cmpd_data)
import numpy as np
new_data = np.array(new_data)
from sklearn.manifold import TSNE
# from sklearn.preprocessing import normalize
# new_data = normalize(new_data,axis=0)
tsne = TSNE( n_components=2, init="pca",  n_iter=2000)
tsne_result = tsne.fit_transform(new_data)
X = tsne_result[:, 0]
Y = tsne_result[:, 1]
for id,compound in enumerate(uniq_compound_list):
    print(compound, X[id],Y[id])
