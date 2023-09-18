#!/usr/bin/env python

import psycopg2
import psycopg2.extras
import traceback
import sys
from openeye.oechem import *

class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, molString, format):
        ifs = oemolistream()
        ifs.SetFormat(format)
        ifs.openstring(molString)
        mol = OEGraphMol()
        OEReadMolecule(ifs,mol)
        return mol

    def convertToInchiKey(self,mol):
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_INCHIKEY)
        ofs.openstring()
        OEWriteMolecule(ofs,mol)
        return ofs.GetString()


class Target:
    def __init__(self,tid, pref_name, target_type, tax_id, organism):
        self.tid = tid
        self.pref_name = pref_name
        self.target_type = target_type
        self.tax_id = tax_id
        self.organism = organism
        self.chemicals = []
        self.components = []

    def addComponent(self, component):
        if component not in self.components:
            self.components.append(component)

    def addChemicals(self,chemical):
        if chemical not in self.chemicals:
            self.chemicals.append(chemical)

    def __eq__(self, other):
        return self.tid == other.tid

class Assay:
    def __init__(self, assay_id, assay_type, description, assay_organism, assay_tax_id,
                 tid, assay_chembl_id):
        self.assay_id = assay_id
        self.assay_type = assay_type
        self.description = description
        self.assay_organism = assay_organism
        self.assay_tax_id = assay_tax_id
        self.tid = tid
        self.assay_chembl_id = assay_chembl_id

    def __eq__(self, other):
        return self.assay_id == other.assay_id

class Component:
    def __init__(self, component_id,uniprot,synonyms,protein_class):
        self.component_id = component_id
        self.uniprot = uniprot
        self.synonyms = synonyms
        self.protein_class = protein_class

    def __eq__(self, other):
        return self.component_id == other.component_id

class Chemical:
    def __init__(self, molregno, smiles, chembl_id):
        self.molregno = molregno
        self.smiles = smiles
        self.chembl_id = chembl_id
        self.activities = []

    def addActivity(self,activity):
        if activity not in self.activities:
            self.activities.append(activity)

    def __eq__(self, other):
        return self.molregno == other.molregno


class Activity:
    def __init__(self,molregno,activity_id,assay_id,standard_type,standard_relation, standard_value, standard_units):
        self.molregno = molregno
        self.activity_id = activity_id
        self.assay_id = assay_id
        self.standard_type = standard_type
        self.standard_relation = standard_relation
        self.standard_value = standard_value
        self.standard_units = standard_units

    def __eq__(self, other):
        return self.activity_id == other.activity_id

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage:%s input.sdf > output.txt"
        sys.exit(1)
    conn = None
    cursor = None
    chemicalsDB = {}
    assayDB = {}
    targetDB = {}
    targetList = []
    try:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        conn = psycopg2.connect(database='chembl_20',user='chembl',host='javelin',password='chembl')
        cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
        while OEReadMolecule(ifs,mol):
            try:
                smiles = OEMolToSmiles(mol)
                cursor.execute("select a.molregno,a.m, b.chembl_id from rdk.mols a inner join chembl_id_lookup b on a.molregno=b.entity_id and b.entity_type='COMPOUND' and m@>%s and m<@%s",(smiles,smiles))
                result = cursor.fetchone()
                if result is not None:
                    molregno = int(result['molregno'])
                    chemical = Chemical(molregno,result['m'], result['chembl_id'])
                    chemicalsDB[molregno] = chemical
            except:
                traceback.print_exc()
                cursor.execute("rollback")

        for molregno in chemicalsDB.keys():
            try:
                cursor.execute("select activity_id, assay_id, molregno, standard_relation, standard_value, standard_units, standard_type from activities where molregno = %s and standard_units='nM' and standard_value < 25",(molregno,))
                result = cursor.fetchall()
                for r in result:
                    activity = Activity(r['molregno'],r['activity_id'],r['assay_id'],r['standard_type'],r['standard_relation'],r['standard_value'],r['standard_units'])
                    chemicalsDB[molregno].addActivity(activity)
                    cursor.execute("select tid from assays where assay_id = %s and assay_organism = 'Homo sapiens'",(r['assay_id'],))
                    r1 = cursor.fetchone()
                    if r1 is not None:
                        cursor.execute("select target_type, pref_name, tax_id, organism from target_dictionary where tid = %s and organism = 'Homo sapiens'",(r1['tid'],))
                        r2 = cursor.fetchone()
                        if r2 is not None:
                            target = Target(r1['tid'],r2['pref_name'],r2['target_type'],r2['tax_id'],r2['organism'])
                            if target not in targetList:
                                target.addChemicals(chemicalsDB[molregno])
                                if target not in targetList:
                                    targetList.append(target)
                                cursor.execute("select component_id from target_components where tid = %s",(r1['tid'],))
                                r3 = cursor.fetchall()
                                for r4 in r3:
                                    component_id = r4['component_id']

                                    cursor.execute("select accession, description, organism from component_sequences where component_id = %s and component_type = 'PROTEIN'",(component_id,))
                                    r5 = cursor.fetchone()
                                    uniprot = r5['accession']

                                    cursor.execute("select component_synonym from component_synonyms where component_id = %s",(component_id,))
                                    r6 = cursor.fetchall()
                                    synonyms = []
                                    for r7 in r6:
                                        synonyms.append(r7['component_synonym'])
                                    target_synonyms = ";".join(synonyms)

                                    cursor.execute("select a.protein_class_id, a.l1,a.l2,a.l3,a.l4,a.l5 from protein_family_classification a inner join component_class b on a.protein_class_id = b.protein_class_id and b.component_id = %s; ",(component_id,))
                                    r9 = cursor.fetchone()
                                    protein_class = "%s\t%s\t%s\t%s\t%s"%(r9['l1'],r9['l2'],r9['l3'],r9['l4'],r9['l5'])
                                    component = Component(component_id,uniprot,target_synonyms, protein_class)
                                    target.addComponent(component)
            except:
                traceback.print_exc()
                cursor.execute("rollback")
        cursor.close()
        ifs.close()
        conn.close()

        for target in targetList:
            if target.target_type=="CELL-LINE":
                continue
            if len(target.components) == 0:
                print "%s\t%s\tNone\tNone\tNone\tNone\tNone\t"%(target.pref_name, target.target_type),
                for c in target.chemicals:
                    print "https://www.ebi.ac.uk/chembl/compound/inspect/%s\t%s"%(c.chembl_id,c.smiles)
                    break
            else:
                for c in target.components:
                    print "%s\t%s\t%s\t"%(c.protein_class, c.synonyms, c.uniprot),
                    for c in target.chemicals:
                        print "https://www.ebi.ac.uk/chembl/compound/inspect/%s\t%s"%(c.chembl_id,c.smiles)
                        break

    except:
        traceback.print_exc()
    finally:
        if cursor is not None:
            cursor.close()
        if conn is not None:
            conn.close()

