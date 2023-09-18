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
from molvs.fragment import LargestFragmentChooser
from molvs.tautomer import TautomerCanonicalizer
from molvs.charge import Uncharger
import pandas


# limit = -1
if __name__ == "__main__":
    sd_reader = Chem.SDMolSupplier("/home/jfeng/Database/HY-L001P_Bioactive Compound Library Plus_Cleaned.sdf")
    #sd_reader = Chem.SDMolSupplier("/home/jfeng/Database/error.sdf")

    #jsonObj = json.loads(open("/home/jfeng/Downloads/JSON export 2019-03-11","r").read())
    try:
        conn = psycopg2.connect(database='medchem_express',user='medchem',host='10.74.2.128',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists chemicals (chemical_id serial primary key not null, molfile text, orig_molfile text, corporate_id varchar(20), mw real, psa real, clogp real);")
        cursor.execute("create schema rdk;")
        # ['fsp3', 'projects', 'composition', 'id', 'molfile', 'class', 'num_rule_of_5_violations', 'isotope_composition', 'molecular_weight', 'log_d',
        #  'isotope_formula', 'name', 'dot_disconnected_formula', 'num_rotatable_bonds', 'p_k_a_type', 'num_h_bond_donors', 'exact_mass',
        #  'batches', 'log_p', 'modified_at', 'heavy_atom_count', 'cns_mpo_score', 'iupac_name', 'created_at', 'topological_polar_surface_area',
        #  'inchi_key', 'synonyms', 'inchi', 'smiles', 'num_h_bond_acceptors', 'owner', 'cxsmiles', 'log_s', 'p_k_a', 'formula']
        progress = 1
        standardizer = Standardizer()
        canonicalizer = TautomerCanonicalizer(max_tautomers=100)
        fragment_chooser = LargestFragmentChooser()
        uncharger = Uncharger()
        for mol in sd_reader:
            if mol is None or mol.GetNumAtoms() == 0:
                continue
            progress += 1
            corporate_id = mol.GetProp("Catalog_NO")
            print(corporate_id)
            orig_molfile = Chem.MolToMolBlock(mol)
            #mol = Chem.MolFromMolBlock(orig_molfile)
            mol.SetProp("_Name",corporate_id)
            mol = fragment_chooser.choose(mol)
            test_mol = uncharger.uncharge(mol)
            problems = Chem.DetectChemistryProblems(test_mol)
            if len(problems) == 0:
                mol = standardizer.standardize(test_mol)
            else:
                print(Chem.MolToSmiles(test_mol), [str(a) for a in problems])
                continue
            mol = canonicalizer.canonicalize(mol)
            molfile = Chem.MolToMolBlock(mol)
            print(Chem.MolToSmiles(mol))
            cursor.execute(sql.SQL("insert into chemicals (molfile,orig_molfile, corporate_id) values (%s, %s,%s)"),[molfile,orig_molfile,corporate_id])
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
