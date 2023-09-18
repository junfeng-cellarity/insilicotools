#!/usr/bin/env python
import csv,sys

if __name__=="__main__":
    if len(sys.argv) != 2:
        print("Usage:%s input.csv"%sys.argv[0])
    else:
        csv_file = open(sys.argv[1],"r")
        csv_reader = csv.DictReader(csv_file)
        for field in csv_reader.fieldnames:
            print(field)