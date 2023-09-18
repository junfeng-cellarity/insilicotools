#!/usr/bin/env python
# compute_tanimoto.py
import sys
import array
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from progressbar import ProgressBar

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
        progress2 = ProgressBar(maxval=len(targets)).start()
        for id2,target in enumerate(targets):
            progress2.update(id2)
            tanimoto = calc_tanimoto(probe,target)
            if tanimoto > best_tanimoto:
                best_tanimoto = tanimoto
        progress2.finish()
        result_list.append(best_tanimoto)
        progress.update(id1)
    progress.finish()



probe = None
target = None
if len(sys.argv) > 2:
    probe = sys.argv[1]
    target = sys.argv[2]
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

    n_processors = 4
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
    for dis in result_list:
        final_result.append(dis)
    stddev = np.std(np.array(final_result))
    mean = np.mean(np.array(final_result))
    n, bins, patches = plt.hist(np.array(final_result),20, normed=True, facecolor='green', alpha=0.5)
    y = mlab.normpdf(bins, mean, stddev)
    plt.plot(bins, y, 'r--')
    plt.xlabel('Tanimoto')
    plt.ylabel('Probability')
    plt.title('Histogram of Similarity Distribution:')

    # Tweak spacing to prevent clipping of ylabel
    # plt.subplots_adjust(left=0.15)
    plt.show()
# print tanimoto(query, target)

else:
    print "usage:%s query.binaryfp target.binaryfp"%sys.argv[0]


