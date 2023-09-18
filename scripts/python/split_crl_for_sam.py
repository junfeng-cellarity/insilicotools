#!/usr/bin/env python
from openpyxl import *

xlsFile = "/Users/jfeng1/forRab/CRL/CRL vs Project Compounds.xlsx"
wb = load_workbook(xlsFile)
sheet = wb.get_active_sheet()
print sheet.rows
bio_number_list = []
project_list = []
for id,row in enumerate(sheet.rows):
    bio_number = str(row[0].value).strip()
    source = str(row[1].value).strip()
    if source == "CRL":
        bio_number_list.append(bio_number)
    if source == "Project":
        project_list.append(bio_number)

print len(bio_number_list),len(project_list)

outputFile = open("/Users/jfeng1/forRab/CRL/bionumber_not_in_project.txt","w")
for bio_number in project_list:
    if bio_number not in bio_number_list:
        print >> outputFile,bio_number
outputFile.close()



