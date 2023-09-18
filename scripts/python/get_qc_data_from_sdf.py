#!/usr/bin/env python
import os,sys
import shutil
from openeye.oechem import *

qc_data_dir = "/home/jfeng/registration/QC_data"
if __name__ == "__main__":
    sdf = sys.argv[1]
    ifs = oemolistream()
    ifs.open(sys.argv[1])
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        if OEHasSDData(mol,"NMR File"):
            nmrfile = OEGetSDData(mol,"NMR File")
            nmrfile_withdir = os.path.join(qc_data_dir,nmrfile)
            shutil.copyfile(nmrfile_withdir,os.path.join(os.getcwd(),nmrfile))
        if OEHasSDData(mol,"LCMS File"):
            lcmsfile = OEGetSDData(mol,"LCMS File")
            lcmsfile_withdir = os.path.join(qc_data_dir,lcmsfile)
            shutil.copyfile(lcmsfile_withdir,os.path.join(os.getcwd(),lcmsfile))
    ifs.close()