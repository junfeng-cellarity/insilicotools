#!/usr/bin/env python
import os,sys
import numpy as np
from openeye.oechem import *

#smarts = "[NX3;!$(N*=O);!$(N-a)]"
smarts_A="[NX4,NX3;!$(N[*]=O);!$(N~a)]"
smarts_B="[NX3;$(NC(=O)OC(C)(C)C)]"
#smarts = "[N;$(N~a)]"
if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Usage:%s input.sdf output.sdf"%sys.argv[0])
    else:
        sub_A = OESubSearch()
        sub_A.Init(smarts_A)

        sub_B = OESubSearch()
        sub_B.Init(smarts_B)

        ofs = oemolostream()
        ofs.open(sys.argv[2])
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        oemol = OEGraphMol()
        while(OEReadMolecule(ifs,oemol)):
            OEDetermineAromaticRingSystems(oemol)
            OEAssignAromaticFlags(oemol)
            OEPrepareSearch(oemol,sub_A)
            distance = []
            for m_A in sub_A.Match(oemol):
                for atm_A in m_A.GetTargetAtoms():
                    coords_A = np.array(oemol.GetCoords(atm_A))
                    print(coords_A)
                    OEPrepareSearch(oemol, sub_B)
                    for m_B in sub_B.Match(oemol):
                        for atm_B in m_B.GetTargetAtoms():
                            coords_B = np.array(oemol.GetCoords(atm_B))
                            print("B",coords_B)
                            distance.append(np.linalg.norm(coords_B-coords_A))
            if len(distance)>0:
                OESetSDData(oemol,"N-N_Distance","%5.3f"%(min(distance)-2.84))
            else:
                OESetSDData(oemol,"N-N_Distance","999.9")
            OEWriteMolecule(ofs,oemol)
        ifs.close()
        ofs.close()