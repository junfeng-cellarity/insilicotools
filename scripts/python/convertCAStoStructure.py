#!/usr/bin/env python
from openeye.oechem import *
from openeye.oemedchem import *
import sys, os
import psycopg2
import psycopg2.extras
import traceback
import re
import requests
import time
if __name__ == "__main__":
    HEADERS = {'user-agent': ('Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_5)'
                              'AppleWebKit/537.36 (KHTML, like Gecko)'
                              'Chrome/45.0.2454.101 Safari/537.36'),
               'referer': 'https://cactus.nci.nih.gov/chemical/structure/'}
    lines = open("/Users/jfeng1/forBrian/cas.list").read().splitlines()
    smiles_list = []
    file = open("/Users/jfeng1/forBrian/cas.smi","w")
    for line in lines:
        cas = line.strip()
        smiles_line = "https://cactus.nci.nih.gov/chemical/structure/%s/smiles" % cas
        response = requests.get(smiles_line, HEADERS)
        file.write("%s %s\n"%(response.content,cas))
    file.close()
