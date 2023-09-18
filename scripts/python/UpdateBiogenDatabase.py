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

database = "/Users/jfeng1/biogen_addon.sdf"
# database = "test.sdf"
# limit = -1
if __name__ == "__main__":
    OEThrow.SetLevel(OEErrorLevel_Error)
    try:
        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='biogen2',user='medchem',host='javelin.biogen.com',password='medchem')
        cursor = conn.cursor()
        n = 0
        bios = {}

        while OEReadMolecule(ifs,mol):
            # if n == limit and limit > 0:
            #     break
            n += 1
            sys.stderr.write("\r%d molecules processed."%n)
            sys.stderr.flush()
            bio_number = OEGetSDData(mol,"BIO Number").strip()
            projname = OEGetSDData(mol, "Concat Distinct;Project")
            if mol.NumAtoms() > 0 and not bios.has_key(bio_number):
                bios[bio_number] = 1
                mol.SetTitle(bio_number)
                molfile = convertMolToString(mol)
                cursor.execute("insert into chemicals (bio_number,molfile,project) select %s,%s,%s where not exists (select bio_number from chemicals where bio_number = %s)",(bio_number,molfile,projname,bio_number))
                cursor.execute("insert into rdk.mols (bio_number,m) select %s, mol_from_ctab(%s) where not exists (select bio_number from rdk.mols where bio_number = %s)",(bio_number,molfile,bio_number))
                cursor.execute("insert into rdk.fps (bio_number, torsionbv, mfp2, ffp2) select %s, torsionbv_fp(mol_from_ctab(%s)), morganbv_fp(mol_from_ctab(%s)), featmorganbv_fp(mol_from_ctab(%s)) where not exists (select bio_number from rdk.fps where bio_number = %s)",(bio_number,molfile,molfile,molfile,bio_number))
        sys.stderr.write("\n")
        cursor.close()
        conn.commit()

    except:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
