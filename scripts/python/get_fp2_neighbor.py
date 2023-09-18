#!/usr/bin/env python
# compute_tanimoto.py
import sys
import array
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from progressbar import ProgressBar
from openeye.oechem import *

# Create a lookup table from byte value to 'on' bit count
# chr(0) has 0 'on' bits
# chr(1) has 1 'on' bit
# chr(2) has 1 'on' bit
# chr(3) has 2 'on' bits
#  ...
# chr(255) has 8 'on' bits
popcount_in_byte = (
    0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
    3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8,
    )

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def calc_tanimoto(query, target):
    a = b = c = 0
    for query_byte, target_byte in zip(query, target):
        a += popcount_in_byte[query_byte]
        b += popcount_in_byte[target_byte]
        c += popcount_in_byte[query_byte & target_byte]
    # I'll define that if neither fingerprint has any
    # set bits (so a == b == 0) then the Tanimoto is 0.
    # Otherwise the equation will try to compute 0/0
    # and Python will throw an exception.
    if not c:
        return 0
    return float(c)/(a+b-c)

def get_best_tanimoto(probes, targets, result_list):
    progress = ProgressBar(maxval=len(probes)).start()
    for id1,probe in enumerate(probes):
        best_tanimoto = 0
        best_id = -1
        progress2 = ProgressBar(maxval=len(targets)).start()
        for id2,target in enumerate(targets):
            progress2.update(id2)
            tanimoto = calc_tanimoto(probe,target)
            if tanimoto > best_tanimoto:
                best_tanimoto = tanimoto
                best_id = id2
        progress2.finish()
        result_list.append((best_tanimoto,best_id))
        progress.update(id1)
    progress.finish()



probe = None
target = None
if len(sys.argv) > 4:
    probe = sys.argv[1]
    target = sys.argv[2]
    probe_sdf = sys.argv[3]
    target_sdf = sys.argv[4]

    mol = OEGraphMol()
    probe_mols = []
    ifs = oemolistream()
    ifs.open(probe_sdf)
    while OEReadMolecule(ifs,mol):
        probe_mols.append(OEGraphMol(mol))
    ifs.close()

    target_mols = []
    ifs = oemolistream()
    ifs.open(target_sdf)
    while OEReadMolecule(ifs,mol):
        target_mols.append(OEGraphMol(mol))
    ifs.close()

    # Use the first fingerprint as the query
    infile = open(probe, "rb")
    s = infile.read(1024//8)  # 128 bytes in a fingerprint
    query = array.array("b", s)

    # Reopen and compute the Tanimoto against all fingerprints
    # including itself.
    infile = open(probe, "rb")
    probes = []
    while 1:
        s = infile.read(1024//8)
        if not s:
            # End of file
            break

        query = array.array("b", s)
        probes.append(query)
    infile.close()

    infile = open(target,"rb")
    targets = []
    while 1:
        s = infile.read(1024//8)
        if not s:
            break
        target = array.array("b", s)
        targets.append(target)

    n_processors = 1
    trunk_size = len(probes)/n_processors
    trunks = list(chunks(probes,trunk_size))
    # import pprint
    # pprint.pprint(trunks)
    import multiprocessing
    manager = multiprocessing.Manager()
    result_list = manager.list()
    jobs = []
    for idx,trunk in enumerate(trunks):
        p = multiprocessing.Process(target=get_best_tanimoto, args=(trunk,targets,result_list))
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()

    final_result = []
    for probe_id,pair in enumerate(result_list):
        target_tanimoto = pair[0]
        target_id = pair[1]
        probe_mol = probe_mols[probe_id]
        target_mol = target_mols[target_id]
        final_result.append((probe_mol,target_mol,target_tanimoto))

    for probe_mol,target_mol,tanimoto in final_result:
        print OEMolToSmiles(probe_mol),OEMolToSmiles(target_mol),tanimoto


else:
    print "usage:%s query.binaryfp target.binaryfp query.sdf target.sdf"%sys.argv[0]


