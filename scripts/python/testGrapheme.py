__author__ = 'jfeng1'
from openeye.oegrapheme import *
from openeye.oedepict import *
from openeye.oechem import *
ifs = oemolistream()
ifs.open("/Users/jfeng1/tmp/template.sdf")
refmol3D = OEGraphMol()
OEReadMolecule(ifs,refmol3D)
ifs.close()

ifs.open("/Users/jfeng1/tmp/target2.sdf")
fitmol3D = OEGraphMol()
OEReadMolecule(ifs,fitmol3D)
ifs.close()

print refmol3D.NumAtoms(),fitmol3D.NumAtoms()

refmol2D = OEGraphMol(refmol3D)
OEPrepareDepictionFrom3D(refmol2D)

fitmol2D = OEGraphMol(fitmol3D)
OEPrepareAlignedDepictionFrom3D(fitmol2D, fitmol3D, refmol2D, refmol3D)


opts = OE2DMolDisplayOptions(300, 300, OEScale_AutoScale)
refscale = OEGetMoleculeScale(refmol2D, opts)
fitscale = OEGetMoleculeScale(fitmol2D, opts)
opts.SetScale(min(refscale, fitscale))

refdisp = OE2DMolDisplay(refmol2D, opts)
OERenderMolecule("OEPrepareAlignedDepictionFrom3D-ref-2D.png", refdisp)

fitdisp = OE2DMolDisplay(fitmol2D, opts)
OERenderMolecule("OEPrepareAlignedDepictionFrom3D-fit-2D.png", fitdisp)
