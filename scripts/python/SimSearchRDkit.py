#!/usr/bin/env python
import rdkit
import os, sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from heapq import heapify, heappush, heappushpop

class MaxHeap():
    def __init__(self, top_n):
        self.h = []
        self.length = top_n
        heapify( self.h)

    def add(self, element):
        if len(self.h) < self.length:
            heappush(self.h, element)
        else:
            heappushpop(self.h, element)

    def getTop(self):
        return sorted(self.h, reverse=True)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage:%s probe.sdf database.sdf result.sdf tag")
        sys.exit(1)
    tag = sys.argv[4]
    probes = [mol for mol in Chem.SDMolSupplier(sys.argv[1]) if mol != None]
    database = [mol for mol in Chem.SDMolSupplier(sys.argv[2]) if mol != None]
    result_dict = {}
    for probe_id,probe in enumerate(probes):
        resultHeap = MaxHeap(5)
        result_dict[probe_id] = resultHeap
        probe_fp = Chem.AllChem.GetMorganFingerprint(probe,3)
        for mol_id, mol in enumerate(database):
            fp = Chem.AllChem.GetMorganFingerprint(mol,3)
            similarity = DataStructs.TanimotoSimilarity(probe_fp,fp)
            resultHeap.add((similarity, mol_id))
    writer = Chem.SDWriter(sys.argv[3])
    result_ids = []
    for key in result_dict:
        print(key, result_dict[key].getTop())
        probe_mol = probes[key]
        for similarity,mol_id in result_dict[key].getTop():
            if mol_id not in result_ids:
                result_mol = database[mol_id]
                result_mol.SetProp("NeighborOf",probe_mol.GetProp(tag))
                result_mol.SetProp("Similarity","%f"%similarity)
#               result_mol.SetProp("Priority",probe_mol.GetProp("Priority"))
                writer.write(result_mol)
                result_ids.append(mol_id)
