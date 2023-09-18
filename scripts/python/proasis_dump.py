#!/usr/bin/env python
import cx_Oracle, paramiko
from scp import SCPClient
import os
import json

def createSSHClient(server, port, user, password):
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server, port, user, password)
    return client


local_proasis_dir = "/Users/jfeng1/Datasets/Kinase/proasis_biib_all"

if __name__ == "__main__":
    ssh = createSSHClient("proasis2.biogen.com",port=22, user="jfeng1",password="Welcome12345")
    scp = SCPClient(ssh.get_transport())
    connection = cx_Oracle.connect("proasis","proasis123","edp02-scan.biogen.com:1725/PRIGPDB")
    cursor = connection.cursor()
#     cursor.execute("""select ts.strucid, ts.structype, ts.strucfile, ts.title, ts.datesolved, ts.depositor, tp.project from tstructure ts
# inner join tprojectstruc tp on ts.strucid=tp.strucid
# where ts.title like '%KINASE%' and ts.title like '%INHIBITOR%' """)
    cursor.execute("""select ts.strucid, ts.structype, ts.strucfile, ts.title, ts.datesolved, ts.depositor, tp.project from tstructure ts 
    inner join tprojectstruc tp on ts.strucid=tp.strucid""")
    for id,type,file,title,date,name, project in cursor:
        dict = {}
        dict['id'] = id
        dict['type'] = type
        dict['path'] = file
        dict['title'] = title
        dict['date'] = date
        dict['name'] = name
        dict['project'] = project
        print json.dumps(dict)
        fname = os.path.basename(file)
        fname = os.path.join(local_proasis_dir,fname)
        scp.get(file,local_path=fname)
    cursor.close()
    connection.close()









