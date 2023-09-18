#!/usr/bin/env python
__author__ = 'jfeng1'
import glob
import os

jarFiles = glob.glob("*.jar")
for jar in jarFiles:
    os.system("jarsigner -keystore /Users/jfeng1/JavaProjects/medchem %s insilico"%jar)
