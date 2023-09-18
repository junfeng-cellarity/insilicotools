#!/usr/bin/env python
import os,sys,subprocess,tempfile
from openeye.oechem import *
OE_OMEGA_COMMAND_PARAM = "%s -mpi_np 8 -in %s -out %s -sdEnergy -warts yes -fromCT -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -sampleHydrogens true -ewindow %f -rms %f "
#OE_OMEGA_COMMAND_PARAM = "%s -mpi_np 8 -in %s -out %s -warts no -fromCT -maxConfRange 200,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600 -rangeIncrement 1 -param %s -addtorlib %s -sampleHydrogens true -ewindow %f -rms %f "
OE_PARAM_DIR = "/home/jfeng/Models/OEModels"
OE_DIR = "/home/jfeng/openeye/bin/"
TAG = "mmff94smod_NoEstat"

def generate_relative_energies_oeomega(input_filename, output_filename):
    ifs = oemolistream()
    ifs.open(input_filename)

    ofs = oemolostream()
    ofs.open(output_filename)

    mcMols = ifs.GetOEMols()
    for mcmol in mcMols:
        lowest_energy = 0
        for conf_id,conf in enumerate(mcmol.GetConfs()):
            energy = float(OEGetSDData(conf, TAG))
            if conf_id == 0:
                lowest_energy = energy
            else:
                if lowest_energy > energy:
                    lowest_energy = energy
        print (lowest_energy)
        for conf in mcmol.GetConfs():
            energy = float(OEGetSDData(conf, TAG))
            print(energy , lowest_energy)
            OESetSDData(conf, "deltaEnergy", "%5.2f" % (energy - lowest_energy))
        OEWriteMolecule(ofs,mcmol)
    ofs.close()
    return


if __name__=="__main__":
    if len(sys.argv)!= 5:
        print("%s input.sdf output.oeb ewindow rmsd"%sys.argv[0])
    else:
        input_fname = sys.argv[1]
        output_fname = sys.argv[2]
        ewindow = float(sys.argv[3])
        rmsd = float(sys.argv[4])
        #param_file = os.path.join(OE_PARAM_DIR,"confgen.param")
        #torsion_file = os.path.join(OE_PARAM_DIR,"my_torlib.txt")
        fd,tmp_file = tempfile.mkstemp(suffix=".oeb",prefix="omega_")
        omega_command = OE_OMEGA_COMMAND_PARAM%(os.path.join(OE_DIR,"omega2"),input_fname,tmp_file, ewindow, rmsd)
        print(omega_command)
        p = subprocess.Popen(omega_command.split(),stdout=subprocess.PIPE)
        p.communicate()
        generate_relative_energies_oeomega(tmp_file,output_fname)







