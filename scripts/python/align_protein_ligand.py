#!/usr/bin/env python
from openeye.oechem import *
import sys,os

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: %s template.pdb target.pdb target_ligand.sdf"
    else:
        targetReceptor = OEGraphMol()
        targetPdb = sys.argv[2]
        ifs = oemolistream()
        ifs.open(targetPdb)
        OEReadMolecule(ifs,targetReceptor)
        ifs.close()

        template = OEGraphMol()
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        OEReadMolecule(ifs,template)
        ifs.close()

        rot = OEDoubleArray(9)
        trans = OEDoubleArray(3)
        alignment = OEGetAlignment(template, targetReceptor)
        rmsd = OERMSD(template,targetReceptor,alignment,True,True, rot, trans)
        OERotate(targetReceptor,rot)
        OETranslate(targetReceptor,trans)

        ofs = oemolostream()
        filename = sys.argv[2]
        filename = os.path.splitext(os.path.basename(filename))[0]
        ofs.open("%s_aligned.pdb" % (filename))
        OEWriteMolecule(ofs,targetReceptor)
        ofs.close()

        ofs = oemolostream()
        filename = sys.argv[3]
        filename = os.path.splitext(os.path.basename(filename))[0]
        ofs.open("%s_aligned.sdf" % (filename))
        targetLigand = OEGraphMol()

        ifs = oemolistream()
        ifs.open(sys.argv[3])
        while OEReadMolecule(ifs,targetLigand):
            OERotate(targetLigand,rot)
            OETranslate(targetLigand,trans)
            OEWriteMolecule(ofs,targetLigand)
        ifs.close()
        ofs.close()

        print rmsd

