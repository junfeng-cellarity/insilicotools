#!/usr/bin/env python

import psycopg2
import traceback
import os,sys
from socket import *
from openeye.oechem import *
import requests
from bs4 import BeautifulSoup

def convertMolToString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_MDL)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    return ofs.GetString()

# limit = -1
if __name__ == "__main__":
    default_directory = None
    if gethostname() == "javelin" or gethostname() == "modeling":
        default_directory = "/DbUCdd/databases/VENDOR_LIBRARIES/"
    else:
        default_directory = "/Users/jfeng1/BiogenDB/Syncom/"


    try:
        database = os.path.join(default_directory,"syncom.sdf")
        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='syncom',user='medchem',host='javelin.biogen.com',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists syncom_chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists syncom_chemicals (chemical_id serial primary key not null, molfile text);")
        cursor.execute("create schema rdk;")
        n = 0
        while OEReadMolecule(ifs,mol):
            # if n == limit and limit > 0:
            #     break
            n += 1
            sys.stderr.write("\r%d molecules processed."%n)
            sys.stderr.flush()
            if mol.NumAtoms() > 0:
                mol.SetTitle("syncom_%d"%n)
                molfile = convertMolToString(mol)
                cursor.execute("insert into syncom_chemicals (molfile) values (%s)",(molfile,))
        sys.stderr.write("\n")
        cursor.execute("select * into rdk.mols from (select chemical_id,mol_from_ctab(molfile::cstring) m from syncom_chemicals) tmp where m is not null;")
        cursor.execute("create index molidx on rdk.mols using gist(m);")
        cursor.execute("alter table rdk.mols add primary key (chemical_id);")
        cursor.execute("select chemical_id,torsionbv_fp(m) as torsionbv, morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;")
        cursor.execute("create index fps_ttbv_idx on rdk.fps using gist(torsionbv);")
        cursor.execute("create index fps_mfp2_idx on rdk.fps using gist(mfp2);")
        cursor.execute("create index fps_ffp2_idx on rdk.fps using gist(ffp2);")
        cursor.execute("alter table rdk.fps add primary key (chemical_id);")
        cursor.close()
        conn.commit()

    except:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
