#!/usr/bin/env python

import requests
import json
import time
import psycopg2
import numpy as np
from numpy.core.multiarray import ndarray
from openeye.oechem import *
import sys


def molFromString(mol_string):
    oemol = OEGraphMol()
    ifs = oemolistream()
    ifs.SetFormat(OEFormat_SDF)
    ifs.openstring(mol_string)
    OEReadMolecule(ifs,oemol)
    ifs.close()
    return oemol

class SkyhawkDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='skyhawk',user='medchem',host='10.74.2.128',password='medchem')
        
    def get_molecule(self,id):
        cursor = None
        try:
            cursor = self.conn.cursor()
            cursor.execute("select orig_molfile,corporate_id from chemicals where chemical_id = %s",(id,))
            row = cursor.fetchone()
            if row is not None:
                return row[0],row[1]
            return None
        finally:
            if cursor is not None:
                cursor.close()
                
    def close(self):
        self.conn.close()
        
    def __del__(self):
        self.conn.close()
        

class ReadoutDef:
    def __init__(self, jsonText):
        dict = json.loads(jsonText)
        self.name = dict['name']
        self.data_type = dict['data_type']
        self.unit_label = dict['nM']
        self.id = dict['id']
    def __str__(self):
        return "%d %s %s %s"%(self.id,self.name,self.data_type,self.unit_label)

kpuu_columns = ["General Protein Binding:Fu,p_mouse","General Protein Binding:Fu,b_mouse","All Full PK:AUClast,p_Mouse_SC","All Full PK:AUClast,b_Mouse_SC"]
sma_desired_columns = ["Protein FL EC50","EC[1.5 Fold]", "CTG CC50", "Splicing FL EC50", "Splicing D7 IC50",
                   "FHOD3 IC50 Off-target", "Selectivity Marker (FHOD3/D7)", "SENP6 IC50 Off-target",
                   "Cytotoxicity marker (SENP6/D7)", "General MDCK-MDR1 - AMRI & Absorption Systems:Efflux", "General MDCK-MDR1 - Pharmaron:Efflux","General MDCK-MDR1 - AMRI & Absorption Systems:Papp A-B",
                   "General MDCK-MDR1 - Pharmaron:Papp A-B", "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
                   "General Nova Ion Channels:hERG IC50",
                   "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
                   "General Solubility:PBS (pH 7.4) (uM)", "General Metabolic Stability - Liver:T1/2_human","General Metabolic Stability - Liver:T1/2_mouse",
                   "Tier 1 Separator",
                   "General Protein Binding:Fu,b_mouse","General Protein Binding:Fu,p_mouse","General Protein Binding:Fu,p_human",
                   "All Full PK:Kp,uu AUC (BALB/c mouse)_Mouse",
                   "All Full PK:CL_Mouse_IV","All Full PK:CL_Mouse_PO","All Full PK:CL_Mouse_SC",
                   "All Full PK:AUCinf,b_Mouse_IV","All Full PK:AUCinf,b_Mouse_PO","All Full PK:AUCinf,b_Mouse_SC",
                   "All Full PK:AUCinf,p_Mouse_IV","All Full PK:AUCinf,p_Mouse_PO","All Full PK:AUCinf,p_Mouse_SC",
                   "All Full PK:AUClast,b_Mouse_IV","All Full PK:AUClast,b_Mouse_PO","All Full PK:AUClast,b_Mouse_SC",
                   "All Full PK:AUClast,p_Mouse_IV","All Full PK:AUClast,p_Mouse_PO","All Full PK:AUClast,p_Mouse_SC",
                   "Tier 2 Separator",
                   "General Nova Ion Channels:Cav1.2 IC50","General Nova Ion Channels:KvLQT1 IC50","General Nova Ion Channels:Late Nav1.5 IC50"
                   ]

htt_desired_columns = ["EC50 mHTT","mHTT EC90", "CC50 CTG1","CC50/EC50","CTG1 CC20","E49-50 EC50","E49b EC50",
                   "FHOD3 EC50 Off-target","PDXDC1 EC50 Off-target","WDR91 EC50 Off-target","SENP6 EC50 Off-target",
                   "General MDCK-MDR1 - AMRI & Absorption Systems:Efflux", "General MDCK-MDR1 - Pharmaron:Efflux","General MDCK-MDR1 - AMRI & Absorption Systems:Papp A-B","General MDCK-MDR1 - Pharmaron:Papp A-B", "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
                   "General Nova Ion Channels:hERG IC50",
                   "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
                   "General Solubility:PBS (pH 7.4) (uM)", "General Metabolic Stability - Liver:T1/2_human","General Metabolic Stability - Liver:T1/2_mouse",
                   "Tier 1 Separator",
                   "General Protein Binding:Fu,b_mouse","General Protein Binding:Fu,p_mouse","General Protein Binding:Fu,p_human",
                   "All Full PK:Kp,uu AUC (BALB/c mouse)_Mouse",
                   "All Full PK:CL_Mouse_IV","All Full PK:CL_Mouse_PO","All Full PK:CL_Mouse_SC",
                   "All Full PK:AUCinf,b_Mouse_IV","All Full PK:AUCinf,b_Mouse_PO","All Full PK:AUCinf,b_Mouse_SC",
                   "All Full PK:AUCinf,p_Mouse_IV","All Full PK:AUCinf,p_Mouse_PO","All Full PK:AUCinf,p_Mouse_SC",
                   "All Full PK:AUClast,b_Mouse_IV","All Full PK:AUClast,b_Mouse_PO","All Full PK:AUClast,b_Mouse_SC",
                   "All Full PK:AUClast,p_Mouse_IV","All Full PK:AUClast,p_Mouse_PO","All Full PK:AUClast,p_Mouse_SC",
                   "Tier 2 Separator",
                   "General Nova Ion Channels:Cav1.2 IC50","General Nova Ion Channels:KvLQT1 IC50","General Nova Ion Channels:Late Nav1.5 IC50"
                   ]


ar_desired_columns = ["- Splicing E1 IC50:22Rv1:Pharmaron","- Splicing E2b/E1 EC50:22Rv1:Pharmaron", "- Splicing E2E3 IC50:22Rv1:Pharmaron", "- Splicing E2b EC50:22Rv1:Pharmaron","- Splicing E2b EC50:22Rv1:Pharmaron", "- Splicing TBP IC50:22Rv1:Pharmaron",
                        "Cytotoxicity Panel-22RV1 (5 days) IC50:22RV1:Pharmaron",
                        "Cytotoxicity Panel-22RV1 (5 days) Max % Inh:22RV1:Pharmaron",
                        "Cytotoxicity Panel-Miapaca-2 (3 days) IC50:Miapaca-2:Pharmaron",
                        "Cytotoxicity Panel-Miapaca-2 (3 days) Max % Inh:Miapaca-2:Pharmaron",
                        "Cytotoxicity Panel-PC-3 (5 days) IC50:PC-3:Pharmaron",
                        "Cytotoxicity Panel-PC-3 (5 days) Max % Inh:PC-3:Pharmaron",
                        "Cytotoxicity Panel-COR-L23 (3 days) IC50:COR-L23:Pharmaron",
                        "Cytotoxicity Panel-COR-L23 (3 days) Max % Inh:COR-L23:Pharmaron",
                        "Cytotoxicity Panel-HT29 (3 days) IC50:HT29:Pharmaron",
                        "Cytotoxicity Panel-HT29 (3 days) Max % Inh:HT29:Pharmaron",
                        "Cytotoxicity Panel-HT29 (5 days) IC50:HT29:Pharmaron",
                        "Cytotoxicity Panel-HT29 (5 days) Max % Inh:HT29:Pharmaron",
                        "Protein 48hr (ICW) IC50 (nM):22RV1:Pharmaron",
                        "Protein 48hr (ICW) Max Inhibition:22RV1:Pharmaron",
                        "Protein 24hr (ICW) IC50 (nM):22RV1:Pharmaron",
                        "Protein 24hr (ICW) Max Inhibition:22RV1:Pharmaron",

                      "General MDCK-MDR1 - AMRI & Absorption Systems:Efflux", "General MDCK-MDR1 - Pharmaron:Efflux",
                      "General MDCK-MDR1 - AMRI & Absorption Systems:Papp A-B", "General MDCK-MDR1 - Pharmaron:Papp A-B",
                      "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
                      "General Nova Ion Channels:hERG IC50",
                      "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
                      "General Solubility:PBS (pH 7.4) (uM)", "General Metabolic Stability - Liver:T1/2_human",
                      "General Metabolic Stability - Liver:T1/2_mouse",
                      "Tier 1 Separator",
                      "General Protein Binding:Fu,b_mouse", "General Protein Binding:Fu,p_mouse",
                      "General Protein Binding:Fu,p_human",
                      "All Full PK:Kp,uu AUC (BALB/c mouse)_Mouse",
                      "All Full PK:CL_Mouse_IV", "All Full PK:CL_Mouse_PO", "All Full PK:CL_Mouse_SC",
                      "All Full PK:AUCinf,b_Mouse_IV", "All Full PK:AUCinf,b_Mouse_PO", "All Full PK:AUCinf,b_Mouse_SC",
                      "All Full PK:AUCinf,p_Mouse_IV", "All Full PK:AUCinf,p_Mouse_PO", "All Full PK:AUCinf,p_Mouse_SC",
                      "All Full PK:AUClast,b_Mouse_IV", "All Full PK:AUClast,b_Mouse_PO",
                      "All Full PK:AUClast,b_Mouse_SC",
                      "All Full PK:AUClast,p_Mouse_IV", "All Full PK:AUClast,p_Mouse_PO",
                      "All Full PK:AUClast,p_Mouse_SC",
                      "Tier 2 Separator",
                      "General Nova Ion Channels:Cav1.2 IC50", "General Nova Ion Channels:KvLQT1 IC50",
                      "General Nova Ion Channels:Late Nav1.5 IC50"
                      ]
'''
'SARS-CoV2 -1PRF', 'SARS-CoV2 anti-viral assay IPK', 
'SARS-CoV2 CTG', 'SARS-CoV2 anti-viral assay single point Evotec', 
'SARS-CoV2-anti-viral-assay-Evotec', 'SARS-CoV2 anti-viral assay Imquest', 
'SARS-CoV2 NMR', 
'SARS-CoV2 anti-viral assay single point Imquest', 'SARS-CoV2 anti-viral assay single point IPK'
'''
covid_desired_columns = ["anti-viral assay single point IPK 10uM SARS-Cov2 average % inhibition:Vero:IPK",
                         "anti-viral assay single point IPK 10 uM SARS-CoV2 average % cell viability:Vero:IPK",
                         "anti-viral assay single point Imquest 10uM SARS-Cov2 average % inhibition:Vero:ImQuest",
                         "anti-viral assay single point Imquest 10 uM SARS-CoV2 average % cell viability:Vero:ImQuest",

                         "-1PRF -1PRF efficiency EC50:HEK293T:Pharmaron",
                         "-1PRF -1PRF efficiency Emax:HEK293T:Pharmaron",
                         "-1PRF Max % Activation:HEK293T:Pharmaron",
                         "-1PRF Max % Inhibition:HEK293T:Pharmaron",

                         "-1PRF -1PRF efficiency EC50:HEK293T:Evotec",
                         "-1PRF -1PRF efficiency Emax:HEK293T:Evotec",
                         "-1PRF Control Emax:HEK293T:Evotec",
                         "-1PRF Control IC50:HEK293T:Evotec",
                         "-1PRF Max % Activation:HEK293T:Evotec",
                         "-1PRF Max % Inhibition:HEK293T:Evotec",
                         "-1PRF SARS Emax:HEK293T:Evotec",
                         "-1PRF SARS IC50:HEK293T:Evotec",
                         "anti-viral assay IPK Calu-3 IC50:Calu-3:IPK",
                         "anti-viral assay IPK Calu-3 CC50:Calu-3:IPK",
                         "anti-viral assay IPK Calu-3 Max (% inhibition):Calu-3:IPK",
                         "anti-viral assay IPK Vero IC50:Vero:IPK",
                         "anti-viral assay IPK Vero CC50:Vero:IPK",
                         "anti-viral assay IPK Vero Max (%inhibition):Vero:IPK",
                         "anti-viral assay Imquest IC50:Vero:Imquest",
                         "anti-viral assay Imquest CC50:Vero:Imquest",
                         "anti-viral assay Imquest Max (%inhibition):Vero:Imquest",
                         "anti-viral assay RetroVirox IC50 w/Pgp inhibitor:Vero:RetroVirox",
                         "anti-viral assay RetroVirox IC90 with Pgp inhibitor:Vero:RetroVirox",
                         "anti-viral assay RetroVirox CC50 w/ Pgp inhibitor:Vero:RetroVirox",
                         "anti-viral assay RetroVirox Conc for max inhibition w/ Pgb inhibitor:Vero:RetroVirox",
                         "anti-viral assay RetroVirox Max(%inhibition) w/ Pgb inhibitor:Vero:RetroVirox"
                         "anti-viral assay RetroVirox IC50:Vero:RetroVirox",
                         "anti-viral assay RetroVirox IC90:Vero:RetroVirox",
                         "anti-viral assay RetroVirox Max (%inhibition):Vero:RetroVirox",
                         "anti-viral assay RetroVirox CC50:Vero:RetroVirox",
                         "anti-viral assay RetroVirox Conc for max inhibition:Vero:RetroVirox",
                         "CTG IC50:GM07491:Pharmaron","CTG IC50:NHDF:Pharmaron",
    "General MDCK-MDR1 - AMRI & Absorption Systems:Efflux", "General MDCK-MDR1 - Pharmaron:Efflux",
    "General MDCK-MDR1 - AMRI & Absorption Systems:Papp A-B", "General MDCK-MDR1 - Pharmaron:Papp A-B",
    "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
    "General Nova Ion Channels:hERG IC50",
    "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
    "General Solubility:PBS (pH 7.4) (uM)", "General Metabolic Stability - Liver:T1/2_human",
    "General Metabolic Stability - Liver:T1/2_mouse",
]

wrn_desired_columns = [
    "- Splicing Splicing E7E8 IC50:A-673:Lake Pharma",
    "- Splicing Splicing E7bE8 EC50:A-673:Lake Pharma",
    "- Splicing Splicing E7E8 IC50:SK-N-MC:Pharmaron",
    "- Splicing Splicing E7bE8 EC50:SK-N-MC:Pharmaron",
    "- Splicing Splicing E14E15 IC50:SK-N-MC:Pharmaron",
    "- Splicing Splicing E13E15 EC50:SK-N-MC:Pharmaron",
    "Cytotoxicity Panel-RKO (3 days) IC50:RKO:Pharmaron",
    "Cytotoxicity Panel-RKO (3 days) Max % Inh:RKO:Pharmaron",
    "Cytotoxicity Panel-RKO (5 days) IC50:RKO:Pharmaron",
    "Cytotoxicity Panel-RKO (5 days) Max % Inh:RKO:Pharmaron",
    "Cytotoxicity Panel-HT29 (3 days) IC50:HT29:Pharmaron",
    "Cytotoxicity Panel-HT29 (3 days) Max % Inh:HT29:Pharmaron",
    "Cytotoxicity Panel-HT29 (5 days) IC50:HT29:Pharmaron",
    "Cytotoxicity Panel-HT29 (5 days) Max % Inh:HT29:Pharmaron",
    "Phototoxicity:IC50 + UV",
    "Phototoxicity:IC50 - UV",
    "Phototoxicity:PIF",
    "General MDCK-MDR1 - AMRI & Absorption Systems:Efflux", "General MDCK-MDR1 - Pharmaron:Efflux",
    "General MDCK-MDR1 - AMRI & Absorption Systems:Papp A-B", "General MDCK-MDR1 - Pharmaron:Papp A-B",
    "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
    "General Nova Ion Channels:hERG IC50",
    "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
    "General Solubility:PBS (pH 7.4) (uM)", "General Metabolic Stability - Liver:T1/2_human",
    "General Metabolic Stability - Liver:T1/2_mouse",
    "Tier 1 Separator",
    "General Protein Binding:Fu,b_mouse", "General Protein Binding:Fu,p_mouse",
    "General Protein Binding:Fu,p_human",
    "All Full PK:Kp,uu AUC (BALB/c mouse)_Mouse",
    "All Full PK:CL_Mouse_IV", "All Full PK:CL_Mouse_PO", "All Full PK:CL_Mouse_SC",
    "All Full PK:AUCinf,b_Mouse_IV", "All Full PK:AUCinf,b_Mouse_PO", "All Full PK:AUCinf,b_Mouse_SC",
    "All Full PK:AUCinf,p_Mouse_IV", "All Full PK:AUCinf,p_Mouse_PO", "All Full PK:AUCinf,p_Mouse_SC",
    "All Full PK:AUClast,b_Mouse_IV", "All Full PK:AUClast,b_Mouse_PO",
    "All Full PK:AUClast,b_Mouse_SC",
    "All Full PK:AUClast,p_Mouse_IV", "All Full PK:AUClast,p_Mouse_PO",
    "All Full PK:AUClast,p_Mouse_SC",
]

cMYB_desired_columns = [
    "- Splicing E11E12 EC50:Jurkat:Pharmaron",
    "- Splicing E11E12 Emax:Jurkat:Pharmaron",
    "- Splicing E11E11b EC50:Jurkat:Pharmaron",
    "- Splicing E11E11b Emax:Jurkat:Pharmaron",
    "- Splicing TBP IC50:Jurkat:Pharmaron",
    "- Splicing TBP Emax:Jurkat:Pharmaron",
    "- Splicing E11E12 EC50:THP-1:Pharmaron",
    "- Splicing E11E11b EC50:THP-1:Pharmaron",
    "- Splicing TBP IC50:THP-1:Pharmaron",
    "- Splicing Selectivity PDXDC1 IC50:Jurkat:Pharmaron",
    "- Splicing Selectivity RAB3IP IC50:Jurkat:Pharmaron",
    "- Splicing Selectivity RNF213 IC50:Jurkat:Pharmaron",
    "- Splicing Selectivity XRN2 IC50:Jurkat:Pharmaron",
    "- Splicing Selectivity RARS1 IC50:Jurkat:Pharmaron",
    "- Splicing Selectivity RBBP8 IC50:Jurkat:Pharmaron",
    "- Splicing Selectivity TBP IC50:Jurkat:Pharmaron",
    "EC50 mHTT",
    "E49b EC50",
    "E49-50 EC50",
    "EC[1.5 Fold]",
    "Protein FL EC50",
    "Splicing D7 IC50",
    "Splicing FL EC50",
    "Cytotoxicity Panel-MOLM13 (3 days) GR50:MOLM13:Pharmaron",
    "Cytotoxicity Panel-MOLM13 (3 days) GRmax:MOLM13:Pharmaron",
    "Cytotoxicity Panel-Jurkat (3 days) GR50:Jurkat:Pharmaron",
    "Cytotoxicity Panel-Jurkat (3 days) GRmax:Jurkat:Pharmaron",
    "Cytotoxicity Panel-AML193 (3 days) GR50:AML193:Pharmaron",
    "Cytotoxicity Panel-AML193 (3 days) GRmax:AML193:Pharmaron",
    "Cytotoxicity Panel-THP1 (3 days) GR50:THP1:Pharmaron",
    "Cytotoxicity Panel-THP1 (3 days) GRmax:THP1:Pharmaron",
    "Cytotoxicity Panel-MV-4-11 (3 days) GR50:MV-4-11:Pharmaron",
    "Cytotoxicity Panel-MV-4-11 (3 days) GRmax:MV-4-11:Pharmaron",
    "Cytotoxicity Panel-HT29 (3 days) GR50:HT29:Pharmaron",
    "Cytotoxicity Panel-HT29 (3 days) GRmax:HT29:Pharmaron",
    "Cytotoxicity Panel-KO52 (3 days) GR50:KO52:Pharmaron",
    "Cytotoxicity Panel-KO52 (3 days) GRmax:KO52:Pharmaron",
    "Cytotoxicity Panel-MEC1 (3 days) GR50:MEC1:Pharmaron",
    "Cytotoxicity Panel-MEC1 (3 days) GRmax:MEC1:Pharmaron",
    "Cytotoxicity Panel-NALM-6 (3 days) GR50:NALM-6:Pharmaron",
    "Cytotoxicity Panel-NALM-6 (3 days) GRmax:NALM-6:Pharmaron",
    "Cytotoxicity Panel-GM07491 (3 days) GR50:GM07491:Pharmaron",
    "Cytotoxicity Panel-GM07491 (3 days) GRmax:GM07491:Pharmaron",
    "Cytotoxicity Panel-NHDF (3 days) GR50:NHDF:Pharmaron",
    "Cytotoxicity Panel-NHDF (3 days) GRmax:NHDF:Pharmaron",
    "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
    "General Nova Ion Channels:hERG IC50",
    "Tier 1 Separator",
    "General Metabolic Stability - Liver:T1/2_human",
    "General Metabolic Stability - Liver:T1/2_mouse",
    "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
    "General Solubility:PBS (pH 7.4) (uM)",
    "Phototoxicity:IC50 + UV",
    "Phototoxicity:IC50 - UV",
    "Phototoxicity:PIF",
    "General MDCK-MDR1 - AMRI & Absorption Systems:Efflux",
    "General MDCK-MDR1 - AMRI & Absorption Systems:Papp A-B",
    "General MDCK-MDR1 - Pharmaron:Efflux",
    "General MDCK-MDR1 - Pharmaron:Papp A-B",
    "General Protein Binding:Fu,b_mouse", "General Protein Binding:Fu,p_mouse",
    "General Protein Binding:Fu,p_human",
    "All Full PK:Kp,uu AUC (BALB/c mouse)_Mouse",
    "All Full PK:CL_Mouse_IV", "All Full PK:CL_Mouse_PO", "All Full PK:CL_Mouse_SC",
    "All Full PK:AUCinf,b_Mouse_IV", "All Full PK:AUCinf,b_Mouse_PO", "All Full PK:AUCinf,b_Mouse_SC",
    "All Full PK:AUCinf,p_Mouse_IV", "All Full PK:AUCinf,p_Mouse_PO", "All Full PK:AUCinf,p_Mouse_SC",
    "All Full PK:AUClast,b_Mouse_IV", "All Full PK:AUClast,b_Mouse_PO",
    "All Full PK:AUClast,b_Mouse_SC",
    "All Full PK:AUClast,p_Mouse_IV", "All Full PK:AUClast,p_Mouse_PO",
    "All Full PK:AUClast,p_Mouse_SC",
]

column_dict = {
    "hERG":["hERG IC50", "hERG IC50 Pharmaron"],
    "mdck":["Efflux", "Papp A-B","Papp B-A"],
    "sma":sma_desired_columns,
    "htt":htt_desired_columns,
    "ar":ar_desired_columns,
    "covid":covid_desired_columns,
    "wrn":wrn_desired_columns,
    "cMYB":cMYB_desired_columns
}


sma_smarts_dict = {
    "W": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1nncs1",
    "Z": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1cncnn1",
    "Y": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1cnccn1",
    "X": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1nnccc1"
}
covid_smarts_dict = {
    "Series 1":"[H][#8]-c1nc[#6,#7;a][#6,#7;a]c1-c1[#6,#7;a][#6,#7;a]c[#6,#7;a][#6,#7;a]1",
    "Series 2":"c1ccc2[#8,#16;a]c(nc2c1)-c1[#6,#7;a][#6,#7;a]c[#6,#7;a][#6,#7;a]1",
    "Series 3":"[#7]S(=O)(=O)c1cccc(c1)-[#6](=O)-[#7]-1-[#6]-[#6]-[#7]-[#6]-[#6]-1"
}

import collections

ar_smarts_dict = {
    "Series 1":"[#6;a]-c1ccc(-[#7]-[#6](=O)-c2cnc(-[#7])cn2)cc1",
    "Series 2":"[#6]-[#7]-[#6](=O)-c1cc(=O)c2cc(-[#7])ccc2o1",
    "Series 3":"[#6;a]-c1cc(=O)n2cc(-[#7])ccc2n1",
    "Series 4":"[#6;a]-[#6,#7;a]1[#6,#7;a][#6,#7;a][#6,#7;a]2[#6,#7;a][#6,#7;a](-[#6,#7;!R0])[#6,#7;a][#6,#7;a]2[#6,#7;a]1"
}

wrn_smarts_dict = {
    "Series 1":"[#6;a]-c1ccc(-[#7]-[#6](=O)-c2cnc(-[#7])cn2)cc1",
    "Series 2":"[#6]-[#7]-[#6](=O)-c1cc(=O)c2cc(-[#7])ccc2o1",
    "Series 3":"[#6;a]-c1cc(=O)n2cc(-[#7])ccc2n1",
    "Series 4":"[#6;a]-[#6,#7;a]1[#6,#7;a][#6,#7;a][#6,#7;a]2[#6,#7;a][#6,#7;a](-[#6,#7;!R0])[#6,#7;a][#6,#7;a]2[#6,#7;a]1"
}

project_smarts_dict = {
    "sma":sma_smarts_dict,
    "htt":sma_smarts_dict,
    "ar":ar_smarts_dict,
    "wrn":wrn_smarts_dict,
    "covid":covid_smarts_dict
}


kpuu_columns = ["General Protein Binding:Fu,p_mouse","General Protein Binding:Fu,b_mouse","All Full PK:AUClast,p_Mouse_SC","All Full PK:AUClast,b_Mouse_SC"]

a = """
- Splicing Splicing E7E8 IC50:A-673:Lake Pharma
- Splicing Splicing E7bE8 EC50:A-673:Lake Pharma
- Splicing Splicing E7cE8 EC50:A-673:Lake Pharma
kpuu
- Splicing Splicing E7bE8 EC50:SK-N-MC:Pharmaron
- Splicing Splicing E7E8 IC50:SK-N-MC:Pharmaron
Cytotoxicity Panel-HT29 Average IC50:HT29:Pharmaron
Cytotoxicity Panel-NCI-H3122 Average IC50:NCI-H3122:Pharmaron
Cytotoxicity Panel-RKO Max % Inh:RKO:Pharmaron
Cytotoxicity Panel-92.1 Average IC50:92.1:Pharmaron
Cytotoxicity Panel-NCI-H3122 IC50:NCI-H3122:Pharmaron
Cytotoxicity Panel-Miapaca-2 Average IC50:Miapaca-2:Pharmaron
Cytotoxicity Panel-HT29 Max % Inh:HT29:Pharmaron
Cytotoxicity Panel-22RV1 Average IC50:22RV1:Pharmaron
Cytotoxicity Panel-92.1 IC50:92.1:Pharmaron
Cytotoxicity Panel-22RV1 IC50:22RV1:Pharmaron
Cytotoxicity Panel-Miapaca-2 IC50:Miapaca-2:Pharmaron
Cytotoxicity Panel-Miapaca-2 Max % Inh:Miapaca-2:Pharmaron
Cytotoxicity Panel-22RV1 Max % Inh:22RV1:Pharmaron
Cytotoxicity Panel-92.1 Max % Inh:92.1:Pharmaron
Cytotoxicity Panel-RKO IC50:RKO:Pharmaron
Cytotoxicity Panel-RKO Average IC50:RKO:Pharmaron
Cytotoxicity Panel-NCI-H3122 Max % Inh:NCI-H3122:Pharmaron
Cytotoxicity Panel-HT29 IC50:HT29:Pharmaron
- Splicing Splicing TBP IC50 (nM):SK-N-MC:Pharmaron
ICW WRN IC50 (nM):SK-N-MC:Pharmaron
ICW WRN Max Inhibition:SK-N-MC:Pharmaron
ICW WRN IC50 (nM):RKO:Pharmaron
ICW WRN Max Inhibition:RKO:Pharmaron
Cytotoxicity Panel-COR-L23 Average IC50:COR-L23:Pharmaron
Cytotoxicity Panel-COR-L23 Max % Inh:COR-L23:Pharmaron
Cytotoxicity Panel-PC-3 Average IC50:PC-3:Pharmaron
Cytotoxicity Panel-COR-L23 IC50:COR-L23:Pharmaron
Cytotoxicity Panel-PC-3 IC50:PC-3:Pharmaron
Cytotoxicity Panel-PC-3 Max % Inh:PC-3:Pharmaron
- Splicing Splicing E3E4 IC50:A-673:Lake Pharma
Cytotoxicity Panel-Lovo Average IC50:Lovo:Pharmaron
Cytotoxicity Panel-Ls1034 IC50:Ls1034:Pharmaron
Cytotoxicity Panel-Lovo IC50:Lovo:Pharmaron
Cytotoxicity Panel-Lovo Max % Inh:Lovo:Pharmaron
Cytotoxicity Panel-Ls1034 Max % Inh:Ls1034:Pharmaron
Cytotoxicity Panel-HCT116 Max % Inh:HCT116:Pharmaron
Cytotoxicity Panel-Ls1034 Average IC50:Ls1034:Pharmaron
Cytotoxicity Panel-HCT116 IC50:HCT116:Pharmaron
Cytotoxicity Panel-HCT116 Average IC50:HCT116:Pharmaron
- Splicing Splicing E7E8 IC50:SK-N-MC:Waltham
- Splicing Splicing E7bE8 EC50:SK-N-MC:Waltham
- Splicing Splicing TBP IC50 (nM):SK-N-MC:Waltham
- Splicing Internal Miniscreen E7E8 Fold Change
- Splicing Internal Miniscreen TBP Fold Change (average)
- Splicing Internal Miniscreen E7E8 Fold Change (average)
- Splicing Internal Miniscreen TBP Fold Change
- Splicing Internal Miniscreen Dose
- Splicing Internal Miniscreen E7bE8 % dCt[IC100] (average)
- Splicing Internal Miniscreen E7bE8 % dCt[IC100]
- Splicing Splicing E7cE8 EC50:HT-29:Lake Pharma
- Splicing SENP6 IC50:HT-29:Lake Pharma
- Splicing PDXDC1 IC50:HT-29:Lake Pharma
- Splicing Splicing E7bE8 EC50:HT-29:Lake Pharma
- Splicing Splicing E7E8 IC50:HT-29:Lake Pharma
- Splicing Splicing E3E4 IC50:HT-29:Lake Pharma
- Cytotoxicity 5 day CTG IC50 (nM):HT-29:Waltham
- Cytotoxicity 5 day CTG IC50 (nM):RKO:Waltham
- Cytotoxicity 5 day CTG IC50 (nM):HT-29:Pharmaron
- Cytotoxicity 3 day CTG IC50 (nM):RKO:Pharmaron
- Cytotoxicity 3 day Max %Inhibition:Lovo:Pharmaron
- Cytotoxicity 5 day CTG IC50 (nM):Ls1034:Pharmaron
- Cytotoxicity 5 day Max %Inhibition:Ls1034:Pharmaron
- Cytotoxicity 3 day CTG IC50 (nM):HT-29:Pharmaron
- Cytotoxicity 3 day Max %Inhibition:RKO:Pharmaron
- Cytotoxicity 3 day Max %Inhibition:Ls1034:Pharmaron
- Cytotoxicity 3 day Max %Inhibition:HT-29:Pharmaron
- Cytotoxicity 3 day CTG IC50 (nM):Lovo:Pharmaron
- Cytotoxicity 3 day CTG IC50 (nM):HCT-116:Pharmaron
- Cytotoxicity 5 day Max %Inhibition:HT-29:Pharmaron
- Cytotoxicity 5 day Max %Inhibition:Lovo:Pharmaron
- Cytotoxicity 5 day CTG IC50 (nM):Lovo:Pharmaron
- Cytotoxicity 5 day CTG IC50 (nM):HCT-116:Pharmaron
- Cytotoxicity 3 day CTG IC50 (nM):Ls1034:Pharmaron
- Cytotoxicity 3 day Max %Inhibition:HCT-116:Pharmaron
- Cytotoxicity 5 day Max %Inhibition:RKO:Pharmaron
- Cytotoxicity 5 day Max %Inhibition:HCT-116:Pharmaron
- Cytotoxicity 5 day CTG IC50 (nM):RKO:Pharmaron
- Splicing eIF Splicing E7E8 IC50 (nM):SK-N-MC:Waltham
- Splicing eIF Splicing E13E15 EC50:SK-N-MC:Waltham
- Splicing eIF Splicing E7cE8 EC50 (nM):SK-N-MC:Waltham
- Splicing Splicing E7cE8 EC50:SK-N-MC:Waltham
- Splicing eIF Splicing E14E15 IC50:SK-N-MC:Waltham
- Splicing Splicing E13E15 EC50:SK-N-MC:Waltham
- Splicing eIF Splicing E7bE8 EC50 (nM):SK-N-MC:Waltham
- Splicing eIF Splicing TBP IC50 (nM):SK-N-MC:Waltham
- Splicing Splicing E14E15 IC50:SK-N-MC:Waltham
- Cytotoxicity 5 day CTG IC50 (nM):HCT-116:Waltham

"""

adme_protocol_names = [
#    "Cytotoxicity Panel",
    "Phototoxicity",
    "General Metabolic Stability - Liver",
    "General hERG",
    "General Solubility",
    "General MDCK-MDR1 - AMRI & Absorption Systems",
    "General MDCK-MDR1 - Pharmaron",
    "General Nova Ion Channels",
    "All Full PK",
    "All Neuro PK",
    "General Protein Binding"
]
       
keywords = ["data_sets","molecules","projects","protocols"]
api_dict = json.loads(requests.request("GET","http://10.74.2.128:8080/data/cdd_api_key.json").text)
cdd_url = "https://app.collaborativedrug.com/api/v1/vaults/4647"
headers = {'X-CDD-token':api_dict['production']}
mol_url = "https://app.collaborativedrug.com/vaults/4647/molecules/%s"

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def is_async_export_finished(id):
    dict = json.loads(requests.request("GET","%s/export_progress/%d"%(cdd_url,id),headers=headers).text)
    if dict['status']=="finished":
        return True
    else:
        return False

def get_cdd_protocols():
    dict = json.loads(requests.request("GET","%s/protocols/?page_size=1000"%(cdd_url), headers=headers).text)
    return dict['objects']

def get_protocol_data(id):
    max_time = 1000
    dict = json.loads(requests.request("GET","%s/protocols/%d/data?async=true"%(cdd_url,id),headers=headers).text)
    export_id = dict['id']
    success = False
    start_time = time.time()
    while True:
        time.sleep(1)
        elapsed_time = time.time()-start_time+1
        if elapsed_time>100:
            elapsed_time = 100
        if elapsed_time>max_time:
            break
        if is_async_export_finished(export_id):
            success = True
            break
    if success:
        return (requests.request("GET","%s/exports/%d"%(cdd_url,export_id),headers=headers)).text
    else:
        return None


protocols = get_cdd_protocols()
project_protocols = []
adme_protocols = []

project_dict = {
    "htt":["HTT Potency"],
    "sma":["SMN2 Potency"],
    "ar":["AR Potency","Cytotoxicity Panel"],
    "hERG":["General hERG"],
    "mdck":["General MDCK-MDR1 - AMRI & Absorption Systems","General MDCK-MDR1 - Pharmaron"],
    "covid":["SARS-CoV2"],
    #"wrn":["WRN Potency","ICW WRN","Cytotoxicity Panel"]
    "wrn": ["WRN Potency","ICW WRN",],
    "cMYB": ["cMYB Potency","HTT Potency","SMN2 Potency","Cytotoxicity Panel"]
}
#"SARS-CoV2 -1PRF","SARS-CoV2 anti-viral assay","SARS-CoV2 anti-viral assay single point","SARS-CoV2 CTG"

if len(sys.argv) != 2:
    print("Usage:%s project_name"%sys.argv[0],file=sys.stderr)
    sys.exit(1)

project_name = sys.argv[1].strip()
if project_name not in project_dict:
    print("Wrong argument: use htt or sma",file=sys.stderr)
    sys.exit(1)


protocol_prefixes = project_dict[project_name]
protocol_names = []
for object in protocols:
    for protocol_prefix in protocol_prefixes:
        if object['name'].startswith(protocol_prefix):
            print("AAA:" + object['name'], protocol_prefix)
            project_protocols.append(object)
            protocol_names.append(object['name'])
            break
    if object['name'] in adme_protocol_names:
        adme_protocols.append(object)
print(protocol_names)
# if True:
#     sys.exit(0)
project_data_dict = {}
project_definitions_dict = {}
for project_protocol in project_protocols:
    protocol_name = project_protocol['name']
    tmp_data = get_protocol_data(project_protocol['id'])
    project_data_dict[protocol_name] = tmp_data
    project_definitions_dict[protocol_name] = project_protocol['readout_definitions']

adme_data_dict = {}
adme_definitions_dict = {}
for adme_protocol in adme_protocols:
    protocol_name = adme_protocol['name']
    tmp_data = get_protocol_data(adme_protocol['id'])
    adme_data_dict[protocol_name] = tmp_data
    adme_definitions_dict[protocol_name] = adme_protocol['readout_definitions']
#print(adme_definitions_dict)
#print(adme_data_dict)

#{'molecule': 53706260, 'readouts': {'413437': {'modifier': '>', 'note': '(n=4)', 'value': 10000.0}, '414205': {'modifier': '<=>', 'value': 1.0}, '413440': {'modifier': '<=>', 'value': 1.0}, '412001': {'modifier': '<=>', 'value': 1.0}, '412003': {'modifier': '>', 'value': 3784.4}, '411999': {'modifier': '>', 'value': 10000.0}, '413438': {'modifier': '>', 'note': '(n=4)', 'value': 10000.0}, '412000': {'modifier': '>', 'value': 10000.0}}, 'batch': 51018110, 'id': 264405693, 'run': 235896}
def_dict = {}
for protocol_name in protocol_names:
    definitions = project_definitions_dict[protocol_name]
    for definition in definitions:
        def_dict[str(definition['id'])] = definition
for adme_name in adme_protocol_names:
    if adme_name in adme_definitions_dict:
        definitions = adme_definitions_dict[adme_name]
        for definition in definitions:
            def_dict[str(definition['id'])] = definition
    else:
        print(adme_definitions_dict.keys())
#print(def_dict)

mol_dict = {}
mol_id_list = []
i = 0
mol_count = 0
mol_data_dict = {}
mol_readout_dict = {}
for protocol_name in project_data_dict.keys():
    data_dict = json.loads(project_data_dict[protocol_name])
    molecules = data_dict['objects']
    for molecule in molecules:
        mol_id = molecule['molecule']
        if mol_id not in mol_id_list:
            mol_id_list.append(mol_id)
            mol_data_dict[mol_id] = {}


for mol_id in mol_id_list:
    for protocol_name in protocol_names:
        protocol_data = json.loads(project_data_dict[protocol_name])
        for mol in protocol_data['objects']:
            if mol['molecule'] == mol_id:
                vendor = None
                vendor_id = -1
                cell_line = None
                cell_line_id = -1
                assay_prefix = None
                if project_name == "covid" or project_name == "wrn" or project_name == "ar" or (project_name == "cMYB" and not (protocol_name.startswith("HTT") or protocol_name.startswith("SMN"))):
                    assay_prefix = protocol_name.replace(project_dict[project_name][0],"").strip()
                    for key in mol['readouts'].keys():
                        if def_dict[key]['name'] == "Vendor":
                            vendor_id = key
                            break
                    if vendor_id in mol['readouts'].keys():
                        vendor = mol['readouts'][vendor_id]['value']
                    for key in mol['readouts'].keys():
                        if def_dict[key]['name'] == "Cell Line":
                            cell_line_id = key
                            break
                    if cell_line_id in mol['readouts'].keys():
                        cell_line = mol['readouts'][cell_line_id]['value']
                for key in mol['readouts'].keys():
                    definition = def_dict[key]
                    def_name = definition['name']
                    def_type = definition['data_type']
                    if def_type == 'Number':
                        if assay_prefix is not None:
                            def_name = "%s %s"%(assay_prefix,def_name)
                        if cell_line is not None:
                            def_name = "%s:%s"%(def_name,cell_line)
                        if vendor is not None:
                            def_name = "%s:%s"%(def_name,vendor)
                        data_point = mol['readouts'][key]
                        if def_name not in mol_data_dict[mol_id]:
                            mol_data_dict[mol_id][def_name] = []
                        if isinstance(data_point,(int,float)):
                            mol_data_dict[mol_id][def_name].append(data_point)
                        elif 'value' in data_point:
                            mol_data_dict[mol_id][def_name].append(data_point['value'])
                        else:
                            continue

    if project_name != "htt" and project_name != "sma" and project_name != "ar" and project_name != "covid" and project_name != "wrn" and project_name != "cMYB":
        continue

    pgp_inhibitor_id = '418460'
    for adme in adme_protocol_names:
        if adme not in adme_data_dict:
            continue
        adme_data = json.loads(adme_data_dict[adme])
        for mol in adme_data['objects']:
            if mol['molecule'] == mol_id:
                specis = None
                specis_id = -1
                for key in mol['readouts'].keys():
                    if def_dict[key]['name'] == "Species":
                        specis_id = key
                        break
                if specis_id in mol['readouts'].keys():
                    specis = mol['readouts'][specis_id]['value']

                RoA = None
                RoA_id = -1
                for key in mol['readouts'].keys():
                    if def_dict[key]['name']=="RoA":
                        RoA_id = key
                        break
                if RoA_id in mol['readouts'].keys():
                    RoA = mol['readouts'][RoA_id]['value']

                pgp_inhibitor = None
                if pgp_inhibitor_id in mol['readouts'].keys():
                    pgp_inhibitor = mol['readouts'][pgp_inhibitor_id]
                for key in mol['readouts'].keys():
                    definition = def_dict[key]
                    def_name = "%s:%s"%(adme,definition['name'])
                    def_type = definition['data_type']
                    if def_type == 'Number':
                        if specis is not None:
                            def_name = "%s_%s" % (def_name, specis)
                        if pgp_inhibitor is not None and pgp_inhibitor == "YES":
                            def_name = "%s_%s" % (def_name, pgp_inhibitor)
                        if RoA is not None:
                            def_name = "%s_%s"%(def_name, RoA)
                        data_point = mol['readouts'][key]
                        if def_name not in mol_data_dict[mol_id]:
                            mol_data_dict[mol_id][def_name] = []
                        if isinstance(data_point,(int,float)):
                            mol_data_dict[mol_id][def_name].append(data_point)
                        elif 'value' in data_point:
                            mol_data_dict[mol_id][def_name].append(data_point['value'])
                        else:
                            continue

molDb = SkyhawkDb()
ofs = oemolostream()
ofs.open("/home/jfeng/apache-tomcat-9.0.13/webapps/data/%s_project.sdf"%project_name)

smarts_dict = None
if project_name in project_smarts_dict:
    smarts_dict = project_smarts_dict[project_name]

for mol_id in mol_data_dict.keys():
    result = molDb.get_molecule(mol_id)
    if result is None:
        print("mol_id %d is not found, please sync the database."%mol_id)
        continue
    molfile,mol_name = result
    oemol = molFromString(molfile)
    OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
    OEAddExplicitHydrogens(oemol)
    if smarts_dict is not None:
        found = False
        for series in smarts_dict.keys():
            ss = OESubSearch(smarts_dict[series])
            OEPrepareSearch(oemol, ss)
            if ss.SingleMatch(oemol):
                OESetSDData(oemol,"Series",series)
                found = True
                break
        if not found:
            OESetSDData(oemol,"Series","Other")
    OESuppressHydrogens(oemol)

    oemol.SetTitle(mol_name)
    OESetSDData(oemol,"mol_id",str(mol_id))
    OESetSDData(oemol,"mol_url",mol_url%str(mol_id))

    if len(column_dict[project_name]) > 0:
        for tag in column_dict[project_name]:
            if tag in mol_data_dict[mol_id].keys():
                    try:
                        values = np.array(mol_data_dict[mol_id][tag])
                        avg_value = np.mean(values)
                        OESetSDData(oemol,tag,"%5.3f"%avg_value)
                    except:
                        OESetSDData(oemol,tag,"")
            else:
                OESetSDData(oemol,tag,"")
    else:
        for tag in mol_data_dict[mol_id].keys():
            try:
                values = np.array(mol_data_dict[mol_id][tag])
                avg_value = np.mean(values)
                OESetSDData(oemol, tag, "%5.3f" % avg_value)
            except:
                OESetSDData(oemol, tag, "")

    tag_list = mol_data_dict[mol_id].keys()
    if kpuu_columns[0] in tag_list and kpuu_columns[1] in tag_list and kpuu_columns[2] in tag_list and kpuu_columns[3] in tag_list:
        try:
            fu_p = np.mean(np.array(mol_data_dict[mol_id][kpuu_columns[0]]))
            fu_b = np.mean(np.array(mol_data_dict[mol_id][kpuu_columns[1]]))
            auc_p = np.mean(np.array(mol_data_dict[mol_id][kpuu_columns[2]]))
            auc_b = np.mean(np.array(mol_data_dict[mol_id][kpuu_columns[3]]))
            OESetSDData(oemol,"kpuu","%5.2f"%(fu_b*auc_b/(fu_p*auc_p)))
        except:
            OESetSDData(oemol,"kpuu","")
    else:
        OESetSDData(oemol, "kpuu", "")
    OEWriteMolecule(ofs,oemol)

ofs.close()



