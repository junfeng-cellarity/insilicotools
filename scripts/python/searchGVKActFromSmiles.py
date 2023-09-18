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

    def getActivityByInchi(self,inchi,name):
        entries = []
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select r.sub_smiles, r.gvk_id, activity_type , activity_value, activity_uom, protein, source, micro_molarvalue from biogen.all_activity_gostar d inner join biogen.structure_details r on d.gvk_id=r.gvk_id where upper(std_activity_type) in ('IC50','Ki','PIC50') and r.inchi = %s;",(inchi,))
            results = cursor.fetchall()
            if results is not None:
                for r in results:
                    try:
                        entry = []
                        entry.append(r['m'])
                        entry.append(name)
                        entry.append(str(r['gvk_id']))
                        entry.append(str(r['activity_type']))
                        entry.append(str(r['activity_uom']).strip())
                        entry.append(str(r['activity_value']))
                        entry.append(str(r['protein']))
                        entry.append(str(r['source']))
                        entry.append(str(r['micro_molarvalue']))
                        entries.append("\t".join(entry))
                    except:
                        traceback.print_exc()
            return entries
        except:
            if cursor is not None:
                cursor.execute("rollback")
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def getActivityBySmiles(self,smiles,name):
        entries = []
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select r.m, r.gvk_id, activity_type , activity_value, activity_uom, protein, source, micro_molarvalue from biogen.all_activity_gostar d inner join rdk.mols r on d.gvk_id=r.gvk_id where upper(std_activity_type) in ('IC50','Ki','PIC50') and m<@mol_from_smiles(%s::cstring) and m@>mol_from_smiles(%s::cstring) and m::text=%s::mol::text;",(smiles,smiles,smiles))
            results = cursor.fetchall()
            if results is not None:
                for r in results:
                    try:
                        entry = []
                        entry.append(r['m'])
                        entry.append(name)
                        entry.append(str(r['gvk_id']))
                        entry.append(str(r['activity_type']))
                        entry.append(str(r['activity_uom']).strip())
                        entry.append(str(r['activity_value']))
                        entry.append(str(r['protein']))
                        entry.append(str(r['source']))
                        entry.append(str(r['micro_molarvalue']))
                        entries.append("\t".join(entry))
                    except:
                        traceback.print_exc()
            return entries
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
    ifs = oemolistream()
    ifs.open("/Users/jfeng1/Datasets/BioActive/BioActiveSet.sdf")
    output = open("/Users/jfeng1/Datasets/BioActive/gvk_out.txt","w")
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        smiles = oechem.OEMolToSmiles(mol)
        print >>sys.stderr,"Searching ",smiles
        results = db.getActivityBySmiles(smiles,mol.GetTitle())
        if results is not None:
            for r in results:
                print >> output, r
    ifs.close()
    output.close()
