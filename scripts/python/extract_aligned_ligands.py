#!/usr/bin/env python
import sys, os
import schrodinger
from schrodinger import structure
from schrodinger.structutils import analyze, transform
import numpy as np

print(schrodinger.get_mmshare_version())
#probe_file = "/Users/jfeng/Projects/Exploratory/probe.maegz"
#aligned_structures_file = "/Users/jfeng/Projects/Exploratory/aligned_kinases.maegz"
if __name__=="__main__":
    if len(sys.argv)!=4:
        print("Usage: %s probe.maegz aligned_structures.maegz aligned_ligands.sdf")
    else:
        probe_file = sys.argv[1]
        probe = structure.StructureReader.read(probe_file)
        aligned_structures = []
        probe_centroid = transform.get_centroid(probe)
        print(probe_centroid)
        aligned_structures_file = sys.argv[2]
        output_ligand_sdf = sys.argv[3]
        with structure.StructureReader(aligned_structures_file) as reader:
            with structure.StructureWriter(output_ligand_sdf) as writer:
                for st in reader:
                    ligands = analyze.find_ligands(st)
                    shortest_distance = 99999
                    closest_st = None
                    for ligand in ligands:
                        lig_centroid = transform.get_centroid(ligand.st)
                        dist = np.linalg.norm(probe_centroid-lig_centroid)
                        if dist < shortest_distance:
                            shortest_distance = dist
                            closest_st = ligand.st
                    if closest_st is not None and shortest_distance < 3.0:
                        writer.append(closest_st)