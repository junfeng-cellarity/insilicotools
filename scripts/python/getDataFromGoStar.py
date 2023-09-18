import glob
import os,sys
import psycopg2
import psycopg2.extras
import traceback
from parse import *
from openeye.oechem import *

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='GoStar',user='medchem',host='javelin',password='medchem')

    def getMolsFromBio(self):
        cursor = None
        molList = []
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select a.*,b.sub_smiles from biogen.all_activity_gostar a inner join biogen.structure_details b on a.gvk_id = b.gvk_id and a.stdname_id in (1107,1108,1109,1110)")
            results = cursor.fetchall()
            for result in results:
                if result is not None:
                    try:
                        keys = result.keys()
                        smiles = result['sub_smiles']
                        mol = OEGraphMol()
                        ifs = oemolistream()
                        ifs.SetFormat(OEFormat_SMI)
                        ifs.openstring(smiles)
                        OEReadMolecule(ifs,mol)
                        mol.SetTitle(str(result['gvk_id']))
                        for key in keys:
                            OESetSDData(mol,key,str(result[key]))
                        molList.append(mol)
                    except:
                        print("Failed to parse %s"%(result['sub_smiles']))
                        continue
            return molList
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    db = BiogenDb()
    molList = db.getMolsFromBio()
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.open("output.sdf")
    for mol in molList:
        if mol is not None:
            OEWriteMolecule(ofs,mol)
    ofs.close()



