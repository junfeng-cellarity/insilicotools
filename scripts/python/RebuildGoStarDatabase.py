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

if __name__ == "__main__":
    try:
        conn = psycopg2.connect(database='GoStar',user='medchem',host='javelin',password='medchem')
        cursor = conn.cursor()
        sys.stderr.write("\n")
        cursor.execute("select * into rdk.mols from (select gvk_id,mol_from_smiles(sub_smiles::cstring) m from biogen.structure_details) tmp where m is not null;")
        cursor.execute("create index molidx on rdk.mols using gist(m);")
        cursor.execute("alter table rdk.mols add primary key (gvk_id);")
        cursor.execute("select gvk_id,torsionbv_fp(m) as torsionbv, morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;")
        cursor.execute("create index fps_ttbv_idx on rdk.fps using gist(torsionbv);")
        cursor.execute("create index fps_mfp2_idx on rdk.fps using gist(mfp2);")
        cursor.execute("create index fps_ffp2_idx on rdk.fps using gist(ffp2);")
        cursor.execute("alter table rdk.fps add primary key (gvk_id);")
        cursor.close()
        conn.commit()

    except:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
