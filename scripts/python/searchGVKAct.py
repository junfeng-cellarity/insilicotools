#!/usr/bin/env python
from openeye.oechem import *
from openeye.oemedchem import *
import sys, os
import psycopg2
import psycopg2.extras
import traceback
import re
class GVKDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='GoStar',user='medchem',host='javelin',password='medchem')

    def getActiveCompouds(self):
        molDict = {}
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select sub_smiles , d.gvk_id, activity_type , activity_value, activity_uom, protein from biogen_201607.all_activity_gostar d inner join biogen_201607.structure_details s on d.gvk_id=s.gvk_id where upper(std_activity_type) in ('IC50','KI','PIC50') and upper(protein) in ('P2X3 RECEPTOR','P2X4 RECEPTOR');")
            results = cursor.fetchall()
            if results is not None:
                for r in results:
                    act_raw_value = r['activity_value']
                    if act_raw_value is None:
                        continue
                    try:
                        smiles = r['sub_smiles']
                        gvk_id = str(r['gvk_id'])
                        act_type = str(r['activity_type']).upper()
                        unit = str(r['activity_uom']).strip()
                        activity = float(act_raw_value)
                        protein = str(r['protein']).upper()
                        # if act_type != 'PIC50':
                        #     if unit != 'uM' and unit != 'nM' and unit!='pM':
                        #         print smiles, activity, unit, " not recorded."
                        #         continue
                        #     if unit == 'uM' and activity > 0.1:
                        #         continue
                        #     if unit == 'nM' and activity > 100:
                        #         continue
                        # else:
                        #     if activity < 6:
                        #         continue
                        if molDict.has_key(gvk_id):
                            mol = molDict[gvk_id]
                        else:
                            mol = OEGraphMol()
                            OEParseSmiles(mol,smiles)
                            mol.SetTitle(gvk_id)
                        OESetSDData(mol,"%s_%s"%(protein,act_type),"%5.2f %s"%(activity,unit))
                        molDict[gvk_id] = mol
                    except:
                        print act_raw_value
                        traceback.print_exc()
                return molDict
        except:
            if cursor is not None:
                cursor.execute("rollback")
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    db = GVKDb()
    molDict = db.getActiveCompouds()
    ofs = oemolostream()
    ofs.open("/Users/jfeng1/gvk_out.sdf")
    for gvk_id in molDict.keys():
        OEWriteMolecule(ofs,molDict[gvk_id])
    ofs.close()
