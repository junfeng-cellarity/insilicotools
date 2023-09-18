#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Usage: %s input.sdf output.sdf"%sys.argv[0])
    else:
        ofs = oemolostream()
        ofs.open(sys.argv[2])

        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        molList = []
        while OEReadMolecule(ifs,mol):
            for dp in OEGetSDDataIter(mol):
                tag = dp.GetTag()
                v = dp.GetValue()
                values = v.strip().split(",")
                if len(values) < 2:
                    continue
                else:
                    for v in values:
                        if len(v) == 0:
                            values.remove(v)
                    new_values = []
                    for v in values:
                        stripped = v.strip()
                        if "requires a valid numeric value.  You entered" in stripped:
                            continue
                        if len(stripped)>0:
                            new_values.append(stripped)
                    dat_value =  ",".join(new_values)
                    OESetSDData(mol,tag,dat_value)
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()