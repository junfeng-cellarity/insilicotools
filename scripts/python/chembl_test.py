__author__ = 'jfeng1'

import MySQLdb
import sys

kinases = ['AAK1',
'BARK',
'CRIK',
'DMPK',
'GRK',
'LATS',
'MAST',
'MRCK',
'NDR',
'PKG',
'PRKY',
'RIOK',
'RSKL',
'SBK1',
'TYRO',
'ULK',
'WNK']

class Target:
    def __init__(self, target_id, target_name):
        self.target_id = target_id
        self.target_name = target_name

    def __eq__(self, other):
        if self.target_id==other.target_id:
            return True
        else:
            return False

con = None
try:
    con = MySQLdb.connect('localhost', 'root', '', 'chembl_20');

    cur = con.cursor()
    targets = []
    for pattern in kinases:
        sql = "select * from target_dictionary where pref_name like binary '%%%s%%'"%pattern
        print sql
        cur.execute(sql)
        result = cur.fetchall()
        print result

except MySQLdb.Error, e:

    print "Error %d: %s" % (e.args[0],e.args[1])
    sys.exit(1)

finally:

    if con:
        con.close()
