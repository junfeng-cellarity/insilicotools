#!/usr/bin/env python
import ftplib,os
FTP_HOST = "ftp-rdb-fr.chem-space.com"
FTP_USER = "guest-user"
FTP_PASS = "BoQKWZdZeDYT"

def save_to_file(ftp,filename):
    local_file = open(filename,"wb")
    ftp.retrbinary("RETR "+filename,local_file.write,1024)
    local_file.close()

directory = "/REAL_Space_June_2020_13.5B_tranches"
ftp = ftplib.FTP(FTP_HOST,FTP_USER,FTP_PASS)
ftp.cwd(directory)

sub_dirs = ['S']
for sub_dir in sub_dirs:
    dir_name = os.path.join(directory, sub_dir)
    ftp.cwd(dir_name)
    dirs = ftp.mlsd()
    for dir in dirs:
        filename = dir[0]
        type = dir[1]['type']
        if type =='dir':
            ftp.cwd(os.path.join(dir_name,filename))
            files = ftp.mlsd()
            for file in files:
                filename = file[0]
                type = file[1]['type']
                if type == "file":
                    save_to_file(ftp,filename)

# ftp.cwd("/REAL_Space_June_2020_13.5B_tranches/S")
# print(ftp.dir())

ftp.quit()