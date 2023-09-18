#!/usr/bin/env python
from openeye.oechem import *
from openeye.oegraphsim import *
import sys, os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import glob
import math
import scipy.stats as stats

OEThrow.SetLevel(OEErrorLevel_Error)
directory = "/Users/junfeng/compound_database/"
properties =   ["ExactMW","SlogP","TPSA (#1)","NumHBD","NumHBA","FractionCSP3"]
propertyIsInt = [False,False,False,True,True,False]
propertyTicks = [[300,350,400,450,500],[-3,-2,-1,0,1,2,3,4,5],[0,25,50,75,100,120,],[0,1,2,3,4,5],[0,1,2,3,4,5,6,7,8,9,10],[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]]
#files = ["wuxi_descriptor.sdf","lifechemical_100k_descriptors.sdf","chemdiv_descriptors.sdf"]
files = ["final_pick.sdf","scaffold_final_sim.sdf"]
#files = ["lifechemical_100k_descriptors.sdf"]
for file in files:
    libraryName = os.path.splitext(file)[0]
    for propertyName in properties:
        plt.clf()
        plt.xlabel(propertyName)
        plt.ylabel('Count')
        plt.title('Histogram of %s Distribution:' % propertyName)
        ifs = oemolistream()
        ifs.open(os.path.join(directory,file))
        mol = OEGraphMol()
        dist = []
        isInt = propertyIsInt[properties.index(propertyName)]
        while OEReadMolecule(ifs,mol):
            sd_data = OEGetSDData(mol, propertyName)
            if(str(sd_data).strip()=="nan" or len(str(sd_data).strip())==0):
                continue
            if isInt:
                property = int(sd_data)
            else:
                property = float(sd_data)
            dist.append(property)

#        print(libraryName, propertyName,len(dist))
        stddev = np.std(np.array(dist))
        mean = np.mean(np.array(dist))
        #weights = np.ones_like(dist) / (len(dist))
        bin = np.array(propertyTicks[(properties.index(propertyName))])
        # print (propertyName,bin)
        # print(mean,stddev)
        # print(bin)
        # print(np.clip(np.array(dist),bin[0],bin[-1]))
        n, bins, patches = plt.hist(np.clip(np.array(dist),bin[0],bin[-1]),bins=bin, density=False, facecolor="blue", alpha=0.4,edgecolor='black', linewidth=1.2)
        print(bins,mean,stddev)
        if not isInt:
            for idx,mybin in enumerate(bins):
                bins[idx] = round(mybin,1)
#        n, bins, patches = plt.hist(np.array(dist),bins=bin, normed=False, facecolor='green', alpha=0.5,edgecolor='black', linewidth=1.2)
        if isInt:
            xtick_labels = (np.array(bins)).astype(int).astype(str)
        else:
            xtick_labels = (np.array(bins)).astype(float).astype(str)
        xtick_labels[0] = "<" + xtick_labels[1]
        xtick_labels[-1] = ">" + xtick_labels[-2]
        plt.xticks(bins, xtick_labels)
        # if propertyName == "hbond_donors" or propertyName == "hbond_acceptors":
        #     xbins = range(int(math.floor(min(dist))), int(math.ceil(max(dist))+1))
        # else:
        #     xbins = bins
        #     bin = len(xbins)
        #     n, bins, patches = plt.hist(np.array(dist),bin, normed=False, facecolor='green', alpha=0.5)
        # else:
        #     n, bins, patches = plt.hist(np.array(dist),10, normed=False, facecolor='green', alpha=0.5)

        #y = stats.norm.pdf(bins, mean, stddev)
        # print("Y=",y)
        # plt.plot(bins, y, 'r--')
        filename = "%s_%s.png" % (libraryName, propertyName)
        print(filename)
        plt.savefig(filename, dpi=300)



    # Tweak spacing to prevent clipping of ylabel
    # plt.subplots_adjust(left=0.15)
    #plt.show()

