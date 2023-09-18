#!/usr/bin/env python
import shlex,sys

datafile = open("/Users/jfeng1/Datasets/Eric/HbF_combined.csv","r")
platefile = open("/Users/jfeng1/Datasets/Eric/control.txt","r")
dict = {}

for lineno,line in enumerate(platefile.readlines()):
    args = line.strip().split('\t')
    # print lineno, args
    if lineno != 0:
        dict[args[0]] = float(args[2])
platefile.close()
#print dict


for lineno,line in enumerate(datafile.readlines()):
    args = line.strip().split(",")
    if lineno == 0:
        #print args
        print "\t".join(args)
    else:
        if len(args)!=15:
            myshlex = shlex.shlex(line.strip())
            myshlex.whitespace = ','
            myshlex.whitespace_split = True
            myshlex.commenters = []
            args = []
            for a in myshlex:
                args.append(a)
            if len(args) != 15:
                print >> sys.stderr,line.strip()
                print >> sys.stderr, args
        raw = float(args[1])
        dmso = float(args[2])
        if raw == dmso:
            plateId = args[6]
            if not plateId.startswith("D"):
                plateId = "D%s"%plateId
            if dict.has_key(plateId):
                dmso = raw/dict[plateId]
                args[2] = "%f"%dmso
                # print >>sys.stderr,"Correctd: before: %f after %f"%(raw,dmso)
            else:
                print >>sys.stderr, "Failed to found:%s,%s"%(args[0],plateId)
        print "\t".join(args)
