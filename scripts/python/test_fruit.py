#!/usr/bin/env python

for a in range(1,1000):
    for b in range(1,1000):
        for c in range(1,1000):
            top = (a*(a+c)*(a+b)+b*(b+c)*(a+b)+c*(a+c)*(b+c))
            bottom = ((a+b)*(a+c)*(b+c))
            if top/bottom == 4 and top%bottom == 0:
                print a,b,c
                break