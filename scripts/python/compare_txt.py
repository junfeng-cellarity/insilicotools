#!/usr/bin/env python

import sys

a = open(sys.argv[1],"r").read().splitlines()
b = open(sys.argv[2],"r").read().splitlines()

duplicate_f = open("duplicate.txt","w")
a_uniq = open("1_uniq.txt","w")
b_uniq = open("2_uniq.txt","w")

a_dict ={}
for a1 in a:
    a_dict[a1] = 1

b_dict = {}
for b1 in b:
    b_dict[b1] = 1

duplicate = []
for a1 in a:
    if b_dict.has_key(a1) and a1 not in duplicate:
        duplicate.append(a1)
for b1 in b:
    if a_dict.has_key(b1) and b1 not in duplicate:
        duplicate.append(b1)

uniq_1 = []
for a1 in a:
    if a1 not in duplicate:
        uniq_1.append(a1)

uniq_2 = []
for b1 in b:
    if b1 not in duplicate:
        uniq_2.append(b1)

print >> a_uniq, "\n".join(uniq_1)
print >> b_uniq, "\n".join(uniq_2)

print >> duplicate_f,"\n".join(duplicate)
