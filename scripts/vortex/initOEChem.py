import sys,glob
import com.dotmatics.vortex.util.Util as Util

os = urllib2 = oechem = None

OE_LICENSE = "openeye/oe_license.txt"

OE_WINDOWS_32_PATH = 'openeye/oejava-2015.Feb.3-Windows-x86.jar'
OE_WINDOWS_64_PATH = 'openeye/oejava-2015.Feb.3-Windows-x64.jar'
OE_MAC_64_PATH = 'openeye/oejava-2015.Feb.3-MacOSX-x64.jar'
OE_LINUX_64_PATH = 'openeye/oejava-2015.Feb.3-Linux-x64.jar'

class VortexInit:

    def __init__(self, vortex):
        self.vortex = vortex
        self.vortexPath = vortex.getVortexFolder()
        self.jarLoader = ClassPathHacker()
        global os
        global urllib2
        import os
        import urllib2

    def appendPathToVortexBase(self, path):
        if path.startswith('/'):
            path = path[1:]
        return os.path.join(self.vortexPath, path)

    def loadJar(self, jar):
        self.jarLoader.addFile(jar)

    def addToPythonPath(self, path, appendToVortexBase=True):
        if appendToVortexBase:
            path = self.appendPathToVortexBase(path)
        if not path in sys.path:
            sys.path.append(path)

    def initJars(self):
        jarLists = open(self.appendPathToVortexBase('jars/jarList.txt'),"r").readlines()
        for jar in jarLists:
            jar = self.appendPathToVortexBase('jars/%s'%jar.strip())
            if jar not in sys.path:
                sys.path.append(jar)
            self.loadJar(jar)
        import chemaxon.license.LicenseManager as LicenseManager
        try:
            LicenseManager.setLicenseFile("http://10.74.2.128:8080/insilico_tools/license.cxl")
        except:
            raise Exception("Failed to find ChemAxon License")


        return

    def initOpeneye(self):
        oechemJar= None
        errorMsg= "Unknown Error"
        if Util.getPlatform() == Util.PlatformIsWindows:
        #    oechemJar = OE_WINDOWS_32_PATH
            if Util.getPlatformArch()==Util.PlatformArchIs64Bit:
                oechemJar = OE_WINDOWS_64_PATH
            elif Util.getPlatformArch()==Util.PlatformArchIs32Bit:
                oechemJar = OE_WINDOWS_32_PATH
            else:
                errorMsg= "Unsupported architecture."
        elif Util.getPlatform() == Util.PlatformIsMac:
            if Util.getPlatformArch()==Util.PlatformArchIs64Bit:
                oechemJar = OE_MAC_64_PATH
            else:
                errorMsg= "Unsupported architecture."
        else:
            oechemJar = OE_LINUX_64_PATH

        if oechemJar:
            print >> sys.stderr, oechemJar
            self.addToPythonPath(oechemJar)
            self.loadJar(oechemJar)
            global oechem
            import openeye.oechem as oechem
            print>>sys.stderr,"OEChem Is Licensed: ",oechem.oechem.OEChemIsLicensed()
            if not oechem.oechem.OEChemIsLicensed():
                licenseString = open(self.appendPathToVortexBase(OE_LICENSE),"r").read();
                oelicense = OELicense(licenseString)
                oelicense.initializeLicense()
                # self.vortex.alert("OEChem Is Licensed: "+ str(oechem.oechem.OEChemIsLicensed()))
        else:
            raise Exception(errorMsg)


    def initRDkit(self):
        #todo:
        #todo: load native libraries using System.load based on operating system.
        #todo:
        return





class OEProduct:
    def __init__(self,productName1):
        self.productName = productName1
        self.licenseKey = None
        self.site = None
        self.licensee = None

    def addLicense(self):
        if self.licenseKey is not None and self.site is not None and self.licensee is not None:
            print >> sys.stderr,"Adding license to %s %d"%(self.productName,oechem.oechem.OEAddLicenseKey(self.licenseKey,self.licensee,self.site))
            return True
        return False

    def toString(self):
        return "%s %s %s %s"%(self.productName,self.licensee,self.site,self.licenseKey)


class OELicense:
    def __init__(self, licenseString):
        self.licenseString = licenseString
        self.oeProducts = []

    def initializeLicense(self):
        for line in self.licenseString.split("\n"):
            if len(line.strip())>0:
                arguments = line.split(":")
                if len(arguments)>=2:
                    key = arguments[0].strip()
                    value = arguments[1].strip()
                    if key.startswith("#PRODUCT"):
                        productName = value
                        oeProduct = OEProduct(productName)
                        self.oeProducts.append(oeProduct)
                        continue
                    elif len(self.oeProducts)==0:
                        continue
                    else:
                        print key,value
                        currentProduct = self.oeProducts[-1]
                        if key.startswith("#LICENSEE"):
                            currentProduct.licensee = value
                        elif key.startswith("#SITE"):
                            currentProduct.site = value
                        elif key.startswith("LICENSE_KEY"):
                            currentProduct.licenseKey = value
                        else:
                            continue
        for product in self.oeProducts:
            product.addLicense()


class ClassPathHacker :
##########################################################
# from http://forum.java.sun.com/thread.jspa?threadID=300557
#
# Author: SG Langer Jan 2007 translated the above Java to this
#       Jython class
# Purpose: Allow runtime additions of new Class/jars either from
#       local files or URL
######################################################
    import java.lang.reflect.Method
    import java.io.File
    import java.net.URL
    import java.net.URLClassLoader
    import jarray

    def addFile (self, s):
        #############################################
        # Purpose: If adding a file/jar call this first
        #       with s = path_to_jar
        #############################################

        # make a URL out of 's'
        f = self.java.io.File (s)
        u = f.toURL ()
        a = self.addURL (u)
        return a

    def addURL (self, u):
        ##################################
        # Purpose: Call this with u= URL for
        #       the new Class/jar to be loaded
        #################################
        sysloader =  self.java.lang.ClassLoader.getSystemClassLoader()
        sysclass = self.java.net.URLClassLoader
        method = sysclass.getDeclaredMethod("addURL", [self.java.net.URL])
        a = method.setAccessible(1)
        jar_a = self.jarray.array([u], self.java.lang.Object)
        b = method.invoke(sysloader, [u])
        return u
