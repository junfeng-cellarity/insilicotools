#!/usr/bin/env python
__author__ = 'jfeng1'

import sys,shlex

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
                cmpdName,dummy,inhibition,geneSymbol,dummy = shlex.split(line.strip())
                targets.append(geneSymbol)
                actDict[geneSymbol] = float(inhibition) #change to 100 - float(inhibition) when needed.
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

        sortedKinaseNames = sorted(standardNameDict.keys())
        text_size = 25
        text_color = "blue"
        for kinase in sortedKinaseNames:
            if newActDict.has_key(kinase):
                act = int("%d"%newActDict[kinase])
                print "at %s"%kinase
                print "scale %d"%act
                print "color red"
                print "circle-lined"
                print "at %s"%kinase
                print "scale %d"%text_size
                print "color %s"%text_color
                print "text %s"%kinase
                print
        print "legend"
        print "color 0 0 0"
        print "next-line"
        print "scale 20"
        print "text 20% inhibition"
        for i in range(9):
            print "space"
        print "color grey"
        print "scale 20"
        print "circle-lined"

        print "scale 40"
        print "next-line"
        print "scale 20"
        print "color black"
        print "text 40% inhibition"
        for i in range(7):
            print "space"
        print "color grey"
        print "scale 40"
        print "circle-lined"

        print "scale 60"
        print "next-line"
        print "scale 20"
        print "color black"
        print "text 60% inhibition"
        for i in range(5):
            print "space"
        print "color grey"
        print "scale 60"
        print "circle-lined"

        print "scale 80"
        print "next-line"
        print "scale 20"
        print "color black"
        print "text 80% inhibition"
        for i in range(3):
            print "space"
        print "color grey"
        print "scale 80"
        print "circle-lined"

        print "scale 100"
        print "next-line"
        print "scale 20"
        print "color black"
        print "text 100% inhibition"
        for i in range(1):
            print "space"
        print "color grey"
        print "scale 100"
        print "circle-lined"

        print "next-line"
        print "color 0 0 0"
        print "scale 20"
        print "text %s"%cmpdName

# print "Names\t%s"%(cmpdName)
        # for kinase in sortedKinaseNames:
        #     act = "-"
        #     if newActDict.has_key(kinase):
        #         act = "%f"%newActDict[kinase]
        #     print "%s\t%s"%(kinase,act)



