#!/usr/bin/env python
from openeye.oechem import *
from openeye.oegraphsim import *
import sys, os
import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

if len(sys.argv)!=3:
    print ("Usage:%s cmpdLibrary.sdf propertyToAnalysis"%sys.argv[0])
else:
    file = sys.argv[1]
    libraryName = os.path.splitext(file)[0]
    propertyName = sys.argv[2]
    ifs = oemolistream()
    ifs.open(file)
    mol = OEGraphMol()
    dist = []
    while OEReadMolecule(ifs,mol):
        property = float(OEGetSDData(mol,propertyName))
        dist.append(property)

    stddev = np.std(np.array(dist))
    mean = np.mean(np.array(dist))
    n, bins, patches = plt.hist(np.array(dist),10, facecolor='green', alpha=0.5)
    y = stats.norm.pdf(bins, mean, stddev)
    plt.plot(bins, y, 'r--')
    plt.xlabel(propertyName)
    plt.ylabel('Count')
    plt.title('Histogram of %s Distribution:'%propertyName)
    plt.savefig("%s_%s.png"%(libraryName,propertyName),dpi=300)

    # Tweak spacing to prevent clipping of ylabel
    # plt.subplots_adjust(left=0.15)
    #plt.show()

