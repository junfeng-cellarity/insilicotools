import pybel
import os,sys

directory = "/Users/jfeng1/Datasets/Arqule"
biogen_fp = os.path.join(directory,"biogen.fpt2")
arqule_fp = os.path.join(directory,"ArQule-FP2.fpt")
pybel.readfile("fpt",arqule_fp)