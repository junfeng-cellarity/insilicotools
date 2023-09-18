import glob
import os,sys
import psycopg2
import psycopg2.extras
import traceback
from parse import *
import json

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='logging',user='medchem',host='javelin',password='medchem')

    def getMessage(self):
        cursor = None
        msgList = []
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select formatted_message from logging_event where caller_method = 'logDocking'")
            results = cursor.fetchall()
            for result in results:
                if result is not None:
                    try:
                        msgList.append(result['formatted_message'])
                    except:
                        print("Failed to parse %s"%(result['formatted_message']))
                        continue
            return msgList
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    db = BiogenDb()
    dict = {}
    msgList = db.getMessage()
    for msg in msgList:
        data = json.loads(msg)
        user =  data['user']
        if not dict.has_key(user):
            dict[user] = 0
        dict[user] = dict[user]+1
    import operator
    sorted_dict = sorted(dict.items(),key=operator.itemgetter(1),reverse=True)
    for v in sorted_dict:
        print str(v[0]),v[1]



