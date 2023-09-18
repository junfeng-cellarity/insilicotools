#!/usr/bin/env python
from openeye.oechem import *
import sys,os

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: %s template.pdb target_1.pdb target_2.pdb ..."
    else:
        numProt = len(sys.argv)-2
        targets = []
        for targetPdb in sys.argv[2:]:
            ifs = oemolistream()
            ifs.open(targetPdb)
            mol = OEGraphMol()
            OEReadMolecule(ifs,mol)
            ifs.close()
            targets.append(mol)

        template = OEGraphMol()
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        OEReadMolecule(ifs,template)
        ifs.close()

        for id,target in enumerate(targets):
            rot = OEDoubleArray(9)
            trans = OEDoubleArray(3)
            alignment = OEGetAlignment(template, target)
            rmsd = OERMSD(template,target,alignment,True,True, rot, trans)
            OERotate(target,rot)
            OETranslate(target,trans)
            ofs = oemolostream()
            filename = sys.argv[id + 2]
            filename = os.path.splitext(os.path.basename(filename))[0]
            ofs.open("%s_out.pdb" % (filename))
            OEWriteMolecule(ofs,target)
            ofs.close()
            print rmsd

