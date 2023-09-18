#!/usr/bin/env python
__author__ = 'jfeng1'
#!/usr/bin/env python
import glob
import os

jarFiles = glob.glob("*.jar")
for jar in jarFiles:
    os.system("jarsigner -tsa http://timestamp.digicert.com -keystore /home/jfeng/Programming/insilicotools/insilicoks -storepass mypassword %s insilico"%jar)
