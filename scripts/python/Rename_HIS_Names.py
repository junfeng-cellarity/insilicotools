#! /usr/bin/env python
# -*- coding: utf8 -*-
""" To rename the residues based on protonization"""
import os,sys
from openeye.oechem import *

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Convert HIS to HID,HIE, HIP based protonation")
        print("Usage:%s input.pdb output.pdb"%sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        protein = OEGraphMol()
        OEReadMolecule(ifs,protein)
        histidineList = []
        for atom in protein.GetAtoms():
            residue = OEAtomGetResidue(atom)
            if residue.GetName() == "HIS":
                if residue not in histidineList:
                    histidineList.append(residue)
        hv = OEHierView(protein)
        for residue in histidineList:
            res = hv.GetResidue(residue.GetChainID(), residue.GetName(), residue.GetResidueNumber())
            #print(residue.GetResidueNumber(), residue.GetName(), residue.GetChainID())
            atom_names = []
            for atom in res.GetAtoms():
                atom_names.append(atom.GetName().strip())
            resname_new="HIS"
            print(atom_names)
            if ("HD1" in atom_names and "HE2" in atom_names):
                resname_new = "HIP"
            elif ("HD1" in atom_names):
                resname_new = "HID"
            elif ("HE2" in atom_names):
                resname_new = "HIE"
            for atom in res.GetAtoms():
                atom_res = OEAtomGetResidue(atom)
                atom_res.SetName(resname_new)
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        OEWriteMolecule(ofs, protein)

# def getatomname(line):
# 	return line[12:16]
#
# def getresname(line):
# 	return line[17:20]
#
# def getresnum(line):
# 	return line[22:26]
#
# if (__name__ == "__main__"):
# 	if (len(sys.argv)!=2):
# 		print "Please give input file(pdb/pqr)."
# 		input()
# 		exit()
# 	fname=sys.argv[1]
# 	fnamelist=os.path.splitext(fname)
# 	fwname=fnamelist[0]+"_HM"+fnamelist[1]
# 	fr=open(fname)
# 	fw=open(fwname,"w")
# 	lines=fr.readlines()
# 	startH=False
# 	endH=0
# 	for i in range(len(lines)):
# 		line=lines[i]
# 		if (line[:4]=="ATOM" or line[:6]=="HETATM"):
# 			# His line have been process
# 			if (startH and i<endH):
# 				continue;
# 			elif (startH and i==endH):
# 				startH=False;
# 			restype=getresname(line);
# 			# first line of HIS residue
# 			if (not startH and (restype=="HIS" or restype=="HID"
# 				or restype=="HIE" or restype=="HIP")):
# 				startH=True
# 				nowrid=getresnum(line)
# 				atomnames=[getatomname(line).strip()]
# 				newlines=[line]
# 				for j in range(i+1,i+30):
# 					line2=lines[j]
# 					rid2=getresnum(line2)
# 					if (rid2!=nowrid):
# 						# set up end pos
# 						endH=j
# 						break
# 					if (not (line2[:4]=="ATOM" or line2[:6]=="HETATM")):
# 						# meet the TER end
# 						endH=j
# 						break
# 					newlines.append(line2)
# 					atomnames.append(getatomname(line2).strip())
# 				# Find the right residue
# 				shouldname="HIS"
# 				if ("HD1" in atomnames and "HE2" in atomnames):
# 					shouldname="HIP"
# 				elif ("HD1" in atomnames):
# 					shouldname="HID"
# 				elif ("HE2" in atomnames):
# 					shouldname="HIE"
# 				for line3 in newlines:
# 					fw.write(line3[:17]+shouldname+line3[20:])
# 				continue
# 		else:
# 			startH=False
# 		fw.write(line)
# 	fr.close()
# 	fw.close()