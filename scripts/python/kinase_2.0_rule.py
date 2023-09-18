#!/usr/bin/env python
__author__ = 'jfeng1'
from openeye.oechem import *
import sys
N_ARO = "[nH0]"
NH_ARO = "[nH1]"
NH_ANILINE = "[c][NH1,NH2]"
CN = "C#N"

def getMatchCount(mIter):
    if mIter is None:
        return 0
    n = 0
    for m in mIter:
        n += 1
    return n

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: %s input.sdf pass.sdf fail.sdf"%sys.argv[0]
    else:
        s_N_ARO = OESubSearch(N_ARO)
        s_NH_ARO = OESubSearch(NH_ARO)
        s_NH_ANILINE = OESubSearch(NH_ANILINE)
        s_CN = OESubSearch(CN)

        passfile = oemolostream()
        passfile.SetFormat(OEFormat_SDF)
        passfile.open(sys.argv[2])

        failfile = oemolostream()
        failfile.SetFormat(OEFormat_SDF)
        failfile.open(sys.argv[3])

        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            OEAddExplicitHydrogens(mol)
            OEPrepareSearch(mol,s_N_ARO)
            m_N_ARO = s_N_ARO.Match(mol,True)
            count_N_ARO = getMatchCount(m_N_ARO)
            OEPrepareSearch(mol,s_NH_ARO)
            count_NH_ARO = getMatchCount(s_NH_ARO.Match(mol,True))
            OEPrepareSearch(mol,s_NH_ANILINE)
            count_NH_ANILINE = getMatchCount(s_NH_ANILINE.Match(mol,True))
            OEPrepareSearch(mol,s_CN)
            count_CN = getMatchCount(s_CN.Match(mol,True))

            if count_N_ARO+count_NH_ARO > 2 or count_CN+count_NH_ANILINE > 0:
                OESuppressHydrogens(mol)
                OEWriteMolecule(passfile,mol)
            else:
                OESuppressHydrogens(mol)
                OEWriteMolecule(failfile,mol)

            # if s.SingleMatch(mol):
            #     failed = True
            #     OESetSDData(mol, "Fail_Reason", name)
            #     OEWriteMolecule(failfile,mol)
            #     break
        ifs.close()
        passfile.close()
        failfile.close()
