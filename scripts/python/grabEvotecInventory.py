#!/usr/bin/env python
from openpyxl import load_workbook
import xlrd
import glob,os
import sqlite3
from sqlitedict import SqliteDict

if __name__ == "__main__":
    inventory_db = SqliteDict("EvotecInventory.db",autocommit=True)
    inventory_db.clear()
    lines = open("/Users/jfeng1/BiogenDB/Evotec/EvotecInventory.csv","r").read().splitlines()
    for lineno,line in enumerate(lines):
        if lineno == 0:
            continue
        args = line.split(",")
        corporate_id = args[1].replace("\"","")
        tmps = corporate_id.rsplit("-", 1)
        if len(tmps)==2:
            bio_number = tmps[0]
            lot_number = tmps[1]
            available_quantity = "%s%s"%(args[4].replace("\"",""),args[5].replace("\"",""))
            value = "%s-%s:%s"%(bio_number,lot_number,available_quantity)
            if inventory_db.has_key(bio_number):
                value = "%s %s"%(inventory_db[bio_number],value)
            inventory_db[bio_number] = value
    print inventory_db["BIO-0878153"]



