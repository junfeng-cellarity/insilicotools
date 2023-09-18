#!/usr/bin/env python
import sys
from openeye.oechem import *
from openeye.oemolprop import *

PSA = "2d PSA"
NROT = "rotatable bonds"
LOGP = "XLogP"
MW = "molecular weight"
hba = "hydrogen-bond acceptors"
hbd = "hydrogen-bond donors"

def hasFluoro(mol):
    for atm in mol.GetAtoms():
        if atm.GetAtomicNum() == 9:
            return True
    return False

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print "Usage:%s input.sdf output.sdf"
    else:
        oechem.OEThrow.SetLevel(oechem.OEErrorLevel_Warning)
        oefilter = OEFilter()
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        ostr = oeosstream()
        oefilter.SetTable(ostr,False)
        propertyNames = ostr.str().strip().split("\t")
        psa_idx = propertyNames.index(PSA)
        mw_idx = propertyNames.index(MW)
        nrot_idx = propertyNames.index(NROT)
        hba_idx = propertyNames.index(hba)
        hbd_idx = propertyNames.index(hbd)
        logp_idx = propertyNames.index(LOGP)
        ostr.clear()
        while OEReadMolecule(ifs,mol):
            if hasFluoro(mol):
                oefilter(mol)
                properties = ostr.str().strip().split("\t")
                psa = float(properties[psa_idx])
                mw = float(properties[mw_idx])
                nrot = float(properties[nrot_idx])
                logp = float(properties[logp_idx])
                hba = float(properties[hba_idx])
                hbd = float(properties[hbd_idx])
                print mw, hbd,hba,logp,nrot,psa
                if mw<300 and hbd<=3 and hba<=3 and logp<=3 and nrot<=3 and psa<=60:
                    print "OK"
                    OEWriteMolecule(ofs,mol)
                ostr.clear()
        ofs.close()





