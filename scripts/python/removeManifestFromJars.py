#!/usr/bin/env python
import glob,os,time,shutil

current_path = os.path.abspath(".")
jars = glob.glob("*.jar")
for jar in jars:
    new_tmp_dir = os.path.join(current_path,"tmp")
    new_dest_dir = os.path.join(current_path,"unsigned")
    new_tmp_jarname = os.path.join(new_tmp_dir,jar)
    new_jarname = os.path.join(new_dest_dir,jar)
    meta_inf = os.path.join(new_tmp_dir,"META-INF")
    shutil.copy(jar,new_tmp_jarname)
    os.chdir(new_tmp_dir)
    os.system("jar xf %s"%new_tmp_jarname)
    shutil.rmtree(meta_inf)
    os.remove(new_tmp_jarname)
    os.system("jar cf %s *"%new_tmp_jarname)
    shutil.copy(new_tmp_jarname,new_jarname)
    tmpfiles = glob.glob(os.path.join(new_tmp_dir,"*"))
    for f in tmpfiles:
        if os.path.isfile(f):
            os.remove(f)
        else:
            if os.path.isdir(f):
                shutil.rmtree(f)
    os.chdir(current_path)
