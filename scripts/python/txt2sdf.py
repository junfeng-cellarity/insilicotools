#!/usr/bin/env python

from openeye.oechem import *
import sys,os
import glob,psycopg2,traceback
import psycopg2.extras
import csv

class SkyhawkDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='skyhawk',user='medchem',host='10.74.2.128',password='medchem')

    def getMolFromDb(self, rgx_number):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select molfile from chemicals where corporate_id = %s",(rgx_number,))
            result = cursor.fetchone()
            if result is not None:
                molfile = result['molfile']
                mol = OEGraphMol()
                ifs = oemolistream()
                ifs.SetFormat(OEFormat_SDF)
                ifs.openstring(molfile)
                OEReadMolecule(ifs,mol)
                return mol

            else:
                return None
        except:
            cursor.execute("rollback")
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()

files = glob.glob("*.csv")
db = SkyhawkDb()
for f in files:
    sdf_fname = "%s.sdf"%os.path.basename(f).split(".")[0]
    ofs = oemolostream()
    ofs.open(sdf_fname)
    titles = []
    title_dict = {}
    csvfile = open(f,"r")
    csv_reader = csv.reader(csvfile)
    n = 0
    for row in csv_reader:
        if n == 0:
            titles = row[:]
            for col_id,title in enumerate(titles):
                title_dict[title] = col_id
            n+=1
            continue
        corporate_id = row[title_dict['ID']]
        mol = db.getMolFromDb(corporate_id)
        mol.SetTitle(corporate_id)
        OESetSDData(mol,"RGX-NUMBER",corporate_id)
        for id,arg in enumerate(row):
            OESetSDData(mol,titles[id],arg)
        OEWriteMolecule(ofs,mol)
    ofs.close()


