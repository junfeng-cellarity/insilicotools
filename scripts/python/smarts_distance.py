#!/usr/bin/env python
import os,sys
import numpy as np
from openeye.oechem import *

#smarts = "[NX3;!$(N*=O);!$(N-a)]"
smarts="[NX4,NX3;!$(N[*]=O);!$(N~a)]"
#smarts = "[N;$(N~a)]"
if __name__=="__main__":
    if len(sys.argv) != 4:
        print("Usage:%s query.sdf input.sdf output.sdf"%sys.argv[0])
    else:
        sub = OESubSearch()
        sub.Init(smarts)

        ifs = oemolistream()
        ifs.open(sys.argv[1])

        qmol = OEGraphMol()
        n_coords = None
        while(OEReadMolecule(ifs,qmol)):
            OEDetermineAromaticRingSystems(qmol)
            OEAssignAromaticFlags(qmol)
            OEPrepareSearch(qmol,sub)
            for m in sub.Match(qmol):
                for atm in m.GetTargetAtoms():
                    n_coords = np.array(qmol.GetCoords(atm))
            break
        ifs.close()
        if n_coords is not None:
            print (n_coords)

        ofs = oemolostream()
        ofs.open(sys.argv[3])                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        ifs.open(sys.argv[2])
        oemol = OEGraphMol()
        while(OEReadMolecule(ifs,oemol)):
            OEDetermineAromaticRingSystems(oemol)
            OEAssignAromaticFlags(oemol)
            OEPrepareSearch(oemol,sub)
            distance = []
            for m in sub.Match(oemol):
                for atm in m.GetTargetAtoms():
                    new_n_coords = np.array(oemol.GetCoords(atm))
                    distance.append(np.linalg.norm(new_n_coords-n_coords))
            if len(distance)>0:
                OESetSDData(oemol,"N-N_Distance","%5.3f"%min(distance))
            else:
                OESetSDData(oemol,"N-N_Distance","999.9")
            OEWriteMolecule(ofs,oemol)


        ifs.close()
        ofs.close()