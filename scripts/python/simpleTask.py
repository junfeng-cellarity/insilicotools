__author__ = 'jfeng1'
#!/usr/bin/env python
import glob
import os

jarFiles = glob.glob("/Users/jfeng1/JavaProjects/vortex/inSilicoTools/lib/*.jar")
for jar in jarFiles:
    jar_name = os.path.basename(jar)
    print "<jar href=\"lib/%s\" />"%jar_name
    os.system("jarsigner -keystore /Users/jfeng1/JavaProjects/medchem %s insilico"%jar)