#!/usr/bin/env python
__author__ = 'jfeng1'

import sys,shlex,parse
from parse import *

TARGETS_DEF = "/Users/jfeng1/KinomeRenderer/TableS1.txt"
if __name__ == "__main__":

    targets = []
    actDict = {}
    cmpdName = None
    standardNameDict = {}
    newActDict = {}

    if len(sys.argv) != 2:
        print "Usage:%s input.txt"%sys.argv[0]
    else:
        inp = sys.argv[1]
        id = 0
        f = open(inp, "r")
        for line in f:
            if id != 0:
                smiles,assayDescription,concentration,unit = line.strip().split("\t")
                geneSymbol = parse("<html><p>PUBCHEM_BIOASSAY: Navigating the Kinome. (Class of assay: other) Panel member name: {}</p></html>",assayDescription)
                print geneSymbol[0],concentration,
                actDict[geneSymbol[0]] = concentration
                targets.append(geneSymbol[0])
            id += 1
        f.close()
        targetDef = open(TARGETS_DEF,"r")
        for lineno,line in enumerate(targetDef.readlines()):
            if lineno == 0:
                continue
            line = line.strip()
            myshlex = shlex.shlex(line)
            myshlex.whitespace = '\t'
            myshlex.whitespace_split = True
            tokens = []
            for a in myshlex:
                tokens.append(a)
            standardNameDict[tokens[1]] = tokens[5]
        targetDef.close()

        for target in targets:
            for key in standardNameDict.keys():
                s = standardNameDict[key]
                if s.lower().find(target.lower()) >= 0:
                    newActDict[key] = actDict[target]
                    break

        print standardNameDict

        for target in targets:
            print target,
            found = False
            for key in standardNameDict.keys():
                s = standardNameDict[key]
                if s == target:
                    print >> sys.stderr, "found %s %s"%(key,target)
                    newActDict[key] = actDict[target]
                    found = True
                    print " found!"
                    break
            if not found:
                print >> sys.stderr, "Not found:",target

        sortedKinaseNames = sorted(standardNameDict.keys())
        print sortedKinaseNames
        print "Names\t%s"%("Kinase Progress Report")
        for kinase in sortedKinaseNames:
            act = "-"
            if newActDict.has_key(kinase):
                act = "%s"%newActDict[kinase]
            print "%s\t%s"%(kinase,act)

