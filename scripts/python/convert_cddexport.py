#!/usr/bin/env python
from openeye.oechem import *
import sys
# if __name__ == "__main__":
#     print ("OK")
# if len(sys.argv)!=3:
#     print ("%s input.sdf output.sdf")
#     print ("Generate a new sdf group by Molecule Name and averaging over all numerical columns")
# else:
def is_float(value):
    new_value = value.replace(">","").replace("<","").replace("=","").strip()
    try:
        float(new_value)
        return True
    except:
        return False

if len(sys.argv)!=3:
    print ("Usage:%s input.sdf output.sdf"%sys.argv[0])
    sys.exit(1)

input = sys.argv[1]
output = sys.argv[2]

ifs = oemolistream()
ifs.open(input)
mol = OEGraphMol()
mol_dict = {}
data_dict = {}
mol_list = []
mol_name_tag = "Molecule Name"
while OEReadMolecule(ifs,mol):
    mol_name = OEGetSDData(mol,mol_name_tag)
    if mol_name not in mol_list:
        mol_list.append(mol_name)
        mol_dict[mol_name] = OEGraphMol(mol)
        data_dict[mol_name] = {}
    for sdtag in OEGetSDDataPairs(mol):
        tag = sdtag.GetTag()
        value = sdtag.GetValue()
        mol_data_dict = data_dict[mol_name]
        if tag not in mol_data_dict:
            mol_data_dict[tag] = []
        if value not in mol_data_dict[tag]:
            mol_data_dict[tag].append(value)

ofs = oemolostream()
ofs.open(output)
for mol_name in mol_list:
    mol = mol_dict[mol_name]
    OEClearSDData(mol)
    for tag in data_dict[mol_name]:
        values = data_dict[mol_name][tag]
        if len(values)>0:
            all_is_float = True
            for v in values:
                if not is_float(v):
                    all_is_float = False
                    break
            if all_is_float:
                sum = 0.0
                for v in values:
                    sum += float(v.replace(">","").replace("<","").replace("=",""))
                sum = sum/len(values)
                OESetSDData(mol,tag,str(sum))
            else:
                OESetSDData(mol,tag,"|".join(values))
        else:
            if len(values)!=0:
                OESetSDData(mol,tag,values[0])
    OEWriteMolecule(ofs,mol)
ofs.close()
