#!/usr/bin/env python

NO_INHIBITORS = "No_inhibitors"
HAS_INHIBITORS = "Has_inhibitors"
NOT_FINISHED = "Not_finished"

# at ADCK3
# scale 10
# color red
# circle-lined

f = open("./data/output.txt","r")
for line in f:
    kinaseName,annotation = line.strip().split()
    if annotation == NO_INHIBITORS:
        print "at %s"%kinaseName
        print "scale 10"
        print "color red"
        print "circle-lined"
        print "text %s"%kinaseName
        print
    if annotation == HAS_INHIBITORS:
        print "at %s"%kinaseName
        print "scale 10"
        print "color green"
        print "circle-lined"
        print "text %s"%kinaseName
        print
    if annotation == NOT_FINISHED:
        print "at %s"%kinaseName
        print "scale 10"
        print "color yellow"
        print "circle-lined"
        print "text %s"%kinaseName
        print
f.close()

