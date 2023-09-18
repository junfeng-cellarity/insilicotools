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
    def __init__(self,molregno,activity_id,assay_id,standard_type,standard_relation, standard_value, standard_units, tid):
        self.molregno = molregno
        self.activity_id = activity_id
        self.assay_id = assay_id
        self.standard_type = standard_type
        self.standard_relation = standard_relation
        self.standard_value = standard_value
        self.standard_units = standard_units
        self.tid = tid

    def __str__(self):
        return "%s_%s_%s"%(self.standard_type.strip(), self.standard_value,self.standard_units.strip())

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
    mols = []
    molDict = {}
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
                newmol = OEGraphMol(mol)
                mols.append(newmol)
                if result is not None:
                    molregno = int(result['molregno'])
                    chemical = Chemical(molregno,result['m'], result['chembl_id'])
                    chemicalsDB[molregno] = chemical
                    molDict[molregno] = newmol
                    OESetSDData(newmol,"CHEMBL", "Found")
                else:
                    OESetSDData(newmol,"CHEMBL","Not Found")
            except:
                traceback.print_exc()
                cursor.execute("rollback")

        for molregno in chemicalsDB.keys():
            try:
                # cursor.execute("select activity_id, assay_id, molregno, standard_relation, standard_value, standard_units, standard_type from activities where molregno = %s and standard_type in ('Ki','Kd')",(molregno,))
                cursor.execute("select a.activity_id, a.assay_id, a.molregno, a.standard_relation, a.standard_value, a.standard_units, a.standard_type from activities a inner join assays b on a.assay_id = b.assay_id and b.assay_type = 'B' and a.standard_type in ('Inhibition','Activity','Residual Activity') and a.molregno = %s",(molregno,))
                result = cursor.fetchall()
                for r in result:
                    value = r['standard_value']
                    if value is None:
                        continue
                    if r['standard_type'] != "Residual Activity":
                        value = 100-value
                    activity = Activity(r['molregno'],r['activity_id'],r['assay_id'],r['standard_type'],r['standard_relation'],value,r['standard_units'], -1)
                    chemicalsDB[molregno].addActivity(activity)
                    # cursor.execute("select tid from assays where assay_id = %s and assay_organism = 'Homo sapiens'",(r['assay_id'],))
                    cursor.execute("select tid from assays where assay_id = %s",(r['assay_id'],))
                    r1 = cursor.fetchone()
                    if r1 is not None:
                        activity.tid = r1['tid']
                        # cursor.execute("select target_type, pref_name, tax_id, organism from target_dictionary where tid = %s and organism = 'Homo sapiens'",(r1['tid'],))
                        cursor.execute("select target_type, pref_name, tax_id, organism from target_dictionary where tid = %s",(r1['tid'],))
                        r2 = cursor.fetchone()
                        if r2 is not None:
                            target = Target(r1['tid'],r2['pref_name'],r2['target_type'],r2['tax_id'],r2['organism'])
                            if target not in targetList:
                                target.addChemicals(chemicalsDB[molregno])
                                targetList.append(target)
                                targetDB[r1['tid']] = target
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
                                        if r7['component_synonym'] is not None:
                                            synonyms.append(r7['component_synonym'])
                                    target_synonyms = ";".join(synonyms)

                                    cursor.execute("select a.protein_class_id, a.l1,a.l2,a.l3,a.l4,a.l5 from protein_family_classification a inner join component_class b on a.protein_class_id = b.protein_class_id and b.component_id = %s; ",(component_id,))
                                    r9 = cursor.fetchone()
                                    protein_class = "%s\t%s\t%s\t%s\t%s"%(r9['l1'],r9['l2'],r9['l3'],r9['l4'],r9['l5'])
                                    component = Component(component_id,uniprot,target_synonyms, protein_class)
                                    target.addComponent(component)
                            else:
                                target = targetDB[r1['tid']]
                                target.addChemicals(chemicalsDB[molregno])
            except:
                traceback.print_exc()
                cursor.execute("rollback")
        cursor.close()
        ifs.close()
        conn.close()

        tagList = []
        tagList.append("CHEMBL")
        for target in targetList:
            if target.target_type=="CELL-LINE":
                continue
            if len(target.components) > 0:
                for c in target.components:
                    synonymList = c.synonyms.split(";")
                    for syn in synonymList:
                        syn = syn.strip().upper()
                        for chemical in target.chemicals:
                            mol = molDict[chemical.molregno]
                            for dp in OEGetSDDataIter(mol):
                                try:
                                    v = float(dp.GetValue())
                                except:
                                    continue
                                if v<30:
                                    if dp.GetTag() not in tagList:
                                        tagList.append(dp.GetTag())
                                oldtag = dp.GetTag().strip().upper()
                                if syn == oldtag:
                                    for act in chemical.activities:
                                        if act.tid == target.tid:
                                            if dp.GetTag() not in tagList:
                                                tagList.append(dp.GetTag())
                                            idx = tagList.index(dp.GetTag())
                                            tag = "%s_Chembl(%%)"%dp.GetTag()
                                            if tag not in tagList:
                                                tagList.insert(idx+1,tag)
                                            actValue = act.standard_value
                                            # if act.standard_type == "uM":
                                            #     actValue = 1000*actValue
                                            OESetSDData(mol, tag, str(actValue))
        print "Smiles\tBIO-NUMBER\t%s"%("\t".join(tagList))
        for mol in mols:
            list = []
            list.append(OEMolToSmiles(mol))
            list.append(OEGetSDData(mol,"Name"))
            for tag in tagList:
                list.append(OEGetSDData(mol,tag)),
            print "\t".join(list)
    except:
        traceback.print_exc()
    finally:
        if cursor is not None:
            cursor.close()
        if conn is not None:
            conn.close()

