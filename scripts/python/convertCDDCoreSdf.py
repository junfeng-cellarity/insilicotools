#!/usr/bin/env python
import os,sys,glob
from openeye.oechem import *

name_dict = {
    "JW":"Jotham Coe",
    "BL":"Brian Lucas",
    "ML":"Michael Luzzio",
    "TW":"Tiansheng Wang"
}
count_dict = {
    "JW":100,
    "BL":100,
    "ML":100,
    "TW":100
}

if __name__ == "__main__":
    sdf_list = glob.glob("*.sdf")
    ofs = oemolostream()
    ofs.open("/home/jfeng/datasets/VirtualIdea/Aug_2019/idea.sdf")
    for sdf in sdf_list:
        key = sdf.split("_", 1)[0]
        chemist = name_dict[key]
        ifs = oemolistream()
        ifs.open(sdf)
        oemol = OEGraphMol()
        while OEReadMolecule(ifs,oemol):
            mol_name = "%s_%d"%(key,count_dict[key])
            count_dict[key] = count_dict[key]+1
            oemol.SetTitle(mol_name)
            OESetSDData(oemol,"Name",mol_name)
            OESetSDData(oemol,"CRO","TBD")
            OESetSDData(oemol,"Chemist",chemist)
            OEWriteMolecule(ofs,oemol)
        ifs.close()
    ofs.close()

    # if len(sys.argv)<2:
    #     sys.exit(1)
    # else:
    #     input_sdf = sys.argv[1]
    #     output_sdf = sys.argv[2]
    #     ofs = oemolostream()
    #     ofs.open(output_sdf)
    #     ifs = oemolistream()
    #     ifs.open(input_sdf)
    #     oemol = OEGraphMol()
    #     while OEReadMolecule(ifs,oemol):
    #         name = OEGetSDData(oemol,"Name")
    #         chemist = ""
    #         if name.startswith("Lib"):
    #             chemist = "Brian Lucas"
    #             name = name.replace("Lib_","Core ")
    #             OESetSDData(oemol,"Name",name)
    #         elif name.startswith("BC"):
    #             chemist = "Brian O'Neil"
    #         elif name.startswith("JF"):
    #             chemist = "Jun Feng"
    #         elif name.startswith("JC"):
    #             chemist = "Jotham Coe"
    #         elif name.startswith("ML"):
    #             chemist = "Michael Luzzio"
    #         else:
    #             chemist = "Unknown"
    #         OESetSDData(oemol,"Chemist",chemist)
    #         OESetSDData(oemol,"Stage","in progress")
    #         OESetSDData(oemol,"CRO","Enamine")
    #         OESetSDData(oemol,"Ranking","Unknown")
    #         oemol.SetTitle(name)
    #         OEWriteMolecule(ofs,oemol)
    #     ifs.close()
    #     ofs.close()
