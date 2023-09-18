#!/usr/bin/env python

import psycopg2
import traceback
import requests
import zipfile
import os
import shutil
import xmlrpc
import json
from rdkit.Chem import Descriptors
from rdkit import Chem
import sys
from psycopg2 import sql
import time
from molvs import Standardizer
from molvs.tautomer import TautomerCanonicalizer

def is_async_export_finished(id):
    dict = json.loads(requests.request("GET","%s/export_progress/%d"%(cdd_url,id),headers=headers).text)
#    print(dict['status'])
    if dict['status']=="finished" or dict['status']=="downloaded":
        return True
    else:
        return False


api_key = "NjE4fHJSTWt2NjQxZ2pHY0ZKSDlkZFRGeUIwc2h4Mm8zTERoZ0xKbXVOZVFSOWJXVVF6Q09RPT0=" #production readonly
cdd_url = "https://app.collaborativedrug.com/api/v1/vaults/5308"
headers = {'X-CDD-token':api_key}

def get_molecules():
    max_time = 100000
    dict = json.loads(requests.request("GET","%s/molecules?async=true"%(cdd_url),headers=headers).text)
    export_id = dict['id']
    success = False
    start_time = time.time()
    while True:
        time.sleep(1)
        elapsed_time = time.time()-start_time
        if elapsed_time>max_time:
            break
        if is_async_export_finished(export_id):
            success = True
            break
    if success:
        jsonTxt = (requests.request("GET","%s/exports/%d"%(cdd_url,export_id),headers=headers)).text
        open("/home/jfeng/virtual_json.txt","w").write(jsonTxt)
        jsonObj = json.loads(jsonTxt)
        return jsonObj

    else:
        return None



# limit = -1
if __name__ == "__main__":
    jsonObj = get_molecules()
    #jsonObj = json.loads(open("/home/jfeng/Downloads/JSON export 2019-03-11","r").read())
    #jsonObj = json.loads(open("/home/jfeng/json.txt","r").read())
    try:
        conn = psycopg2.connect(database='skyhawk_virtual',user='medchem',host='10.74.2.128',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists chemicals (chemical_id integer primary key not null, molfile text, corporate_id varchar(20), mw real, psa real, clogp real);")
        cursor.execute("create schema rdk;")
        # ['fsp3', 'projects', 'composition', 'id', 'molfile', 'class', 'num_rule_of_5_violations', 'isotope_composition', 'molecular_weight', 'log_d',
        #  'isotope_formula', 'name', 'dot_disconnected_formula', 'num_rotatable_bonds', 'p_k_a_type', 'num_h_bond_donors', 'exact_mass',
        #  'batches', 'log_p', 'modified_at', 'heavy_atom_count', 'cns_mpo_score', 'iupac_name', 'created_at', 'topological_polar_surface_area',
        #  'inchi_key', 'synonyms', 'inchi', 'smiles', 'num_h_bond_acceptors', 'owner', 'cxsmiles', 'log_s', 'p_k_a', 'formula']
        progress = 1
        standardizer = Standardizer()
        canonicalizer = TautomerCanonicalizer(max_tautomers=100)
        for molObj in jsonObj['objects']:
            progress+=1
            corporate_id = None
            names = molObj['synonyms']
            chemical_id = molObj['id']
            for name in names:
                if name.startswith("VIR-"):
                    corporate_id = name.strip()
                    break
            if corporate_id is not None:
                molfile = molObj['molfile']
                mol = Chem.MolFromMolBlock(molfile)
                mol = standardizer.standardize(mol)
                mol = canonicalizer.canonicalize(mol)
                molfile = Chem.MolToMolBlock(mol)
                mw = Descriptors.MolWt(mol)
                clogp = Descriptors.rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
                psa = Descriptors.rdMolDescriptors.CalcTPSA(mol)
                cursor.execute(sql.SQL("insert into chemicals (chemical_id, molfile,corporate_id,psa,mw,clogp) values (%s, %s,%s,%s,%s,%s)"),[chemical_id,molfile,corporate_id,psa,mw,clogp])
            else:
                print(molObj)
                print("Failed to find VIR number")
        cursor.execute("select * into rdk.mols from (select chemical_id,mol_from_ctab(molfile::cstring) m from chemicals) tmp where m is not null;")
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
