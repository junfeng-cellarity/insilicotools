#!/usr/bin/env python

import psycopg2
import traceback
import sys
from openeye.oechem import *

def convertMolToString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_MDL)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    return ofs.GetString()

database = "/DbUCdd/databases/VENDOR_LIBRARIES/Whole_internal_database.sdf.gz"
# limit = -1
if __name__ == "__main__":
    try:
        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='biogen',user='medchem',host='javelin',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists chemicals (bio_number char(20) primary key not null, molfile text, project char(50));")
        cursor.execute("create schema rdk;")
        n = 0
        bios = {}
        while OEReadMolecule(ifs,mol):
            # if n == limit and limit > 0:
            #     break
            n += 1
            sys.stderr.write("\r%d molecules processed."%n)
            sys.stderr.flush()
            bio_number = OEGetSDData(mol,"CORPORATE_ID").strip()
            projname = OEGetSDData(mol, "Concat Distinct;Project")
            if mol.NumAtoms() > 0 and not bios.has_key(bio_number):
                bios[bio_number] = 1
                mol.SetTitle(bio_number)
                molfile = convertMolToString(mol)
                cursor.execute("insert into chemicals (bio_number,molfile,project) values (%s,%s,%s)",(bio_number,molfile,projname))
        sys.stderr.write("\n")
        cursor.execute("select * into rdk.mols from (select bio_number,mol_from_ctab(molfile::cstring) m from chemicals) tmp where m is not null;")
        cursor.execute("create index molidx on rdk.mols using gist(m);")
        cursor.execute("alter table rdk.mols add primary key (bio_number);")
        cursor.execute("select bio_number,torsionbv_fp(m) as torsionbv, morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;")
        cursor.execute("create index fps_ttbv_idx on rdk.fps using gist(torsionbv);")
        cursor.execute("create index fps_mfp2_idx on rdk.fps using gist(mfp2);")
        cursor.execute("create index fps_ffp2_idx on rdk.fps using gist(ffp2);")
        cursor.execute("alter table rdk.fps add primary key (bio_number);")
        cursor.close()
        conn.commit()

    except:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
