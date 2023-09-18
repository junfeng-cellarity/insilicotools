#!/usr/bin/env python

import psycopg2
import traceback
import sys
from openeye.oechem import *
from rdkit.Chem import Descriptors
from rdkit import Chem
from psycopg2 import sql

def calculateCLogP_PSA(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        return logp,tpsa
    except:
        return 0,0

def convertMolToString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_MDL)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    return ofs.GetString()

database = "/home/jfeng/Database/AldrichMarketSelect/MarketSelectMain2019week0105.sdf"
# limit = -1
if __name__ == "__main__":
    try:
        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='market_select',user='medchem',host='10.74.2.128',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists aldrich_chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists aldrich_chemicals (structure_id char(20) primary key not null, molfile text, buildingblock char(10), screeningcompound char(10), mw real, psa real, clogp real);")
        cursor.execute("create schema rdk;")
        n = 0
        dict = {}
        while OEReadMolecule(ifs,mol):
            # if n == limit and limit > 0:
            #     break
            n += 1
            # if n > 10000:
            #     break
            sys.stderr.write("\r%d molecules processed."%n)
            sys.stderr.flush()
            structure_id = OEGetSDData(mol,"STRUCTURE_ID").strip()
            bb = OEGetSDData(mol, "BB")
            sc = OEGetSDData(mol, "SC")
            if mol.NumAtoms() > 0 and structure_id not in dict:
                dict[structure_id] = 1
                mol.SetTitle(structure_id)
                molfile = convertMolToString(mol).decode('UTF-8')
                OEDeleteEverythingExceptTheFirstLargestComponent(mol)
                mw = OECalculateMolecularWeight(mol)
                smiles = OEMolToSmiles(mol)
                try:
                    rd_mol = Chem.MolFromSmiles(smiles)
                    clogp = Descriptors.MolLogP(rd_mol)
                    psa = Descriptors.TPSA(rd_mol)
    #               cursor.execute("insert into aldrich_chemicals (structure_id,molfile,buildingblock,screeningcompound,mw,psa,clogp) values (%s,%s,%s,%s,%s,%s,%s)",(structure_id,molfile,bb,sc,mw,psa,clogp))
                    cursor.execute(sql.SQL("insert into aldrich_chemicals (structure_id,molfile,buildingblock,screeningcompound,mw,psa,clogp) values (%s,%s,%s,%s,%s,%s,%s)"),[structure_id,molfile,bb,sc,mw,psa,clogp])
                except:
                    print("Illegal smiles: %s"%smiles,file=sys.stderr)
        sys.stderr.write("\n")
        cursor.execute("select * into rdk.mols from (select structure_id,mol_from_ctab(molfile::cstring) m from aldrich_chemicals) tmp where m is not null;")
        cursor.execute("create index molidx on rdk.mols using gist(m);")
        cursor.execute("alter table rdk.mols add primary key (structure_id);")
        cursor.execute("select structure_id,torsionbv_fp(m) as torsionbv, morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;")
        cursor.execute("create index fps_ttbv_idx on rdk.fps using gist(torsionbv);")
        cursor.execute("create index fps_mfp2_idx on rdk.fps using gist(mfp2);")
        cursor.execute("create index fps_ffp2_idx on rdk.fps using gist(ffp2);")
        cursor.execute("alter table rdk.fps add primary key (structure_id);")
        cursor.close()
        conn.commit()

    except:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
