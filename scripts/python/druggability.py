#!/usr/bin/env python

from openeye.oechem import *
from openeye.oespicoli import *
from openeye.oegrid import *
import numpy
import math,operator,sys

def readInFacetPoints():
    facetPoints= [(0.20358218145599999, 0.19457979239500001, 0.103406961977), (0.15771933363599999, 0.25519524250600001, 0.0), (0.0, 0.15771933363599999, -0.25519524250600001), (-0.15771933363599999, 0.25519524250600001, 0.0), (0.0, 0.15771933363599999, 0.25519524250600001), (0.25519524250600001, 0.0, 0.15771933363599999), (-0.25519524250600001, 0.0, 0.15771933363599999), (0.0, -0.15771933363599999, 0.25519524250600001), (-0.25519524250600001, 0.0, -0.15771933363599999), (-0.15771933363599999, -0.25519524250600001, 0.0), (0.25519524250600001, 0.0, -0.15771933363599999), (0.15771933363599999, -0.25519524250600001, 0.0), (0.0, -0.15771933363599999, -0.25519524250600001), (0.041995361459999998, 0.042804523999799997, -0.29394584933500001), (-0.20358218145599999, 0.19457979239500001, -0.103406961977), (-0.041995361459999998, 0.042804523999799997, 0.29394584933500001), (0.24557754291600001, 0.168625231643, -0.035457039764399997), (-0.078261563764200001, 0.27203219361999997, -0.099366056939399997), (-0.12532061769200001, 0.14621148597600001, 0.23003683216000001), (0.23003683216000001, -0.12532061769200001, 0.14621148597600001), (-0.27203219361999997, -0.099366056939399997, 0.078261563764200001), (0.14621148597600001, 0.23003683216000001, -0.12532061769200001), (-0.035457039764399997, 0.24557754291600001, 0.168625231643), (-0.103406961977, -0.20358218145599999, 0.19457979239500001), (0.29394584933500001, -0.041995361459999998, 0.042804523999799997), (0.103406961977, 0.20358218145599999, 0.19457979239500001), (-0.099366056939399997, -0.078261563764200001, 0.27203219361999997), (0.099366056939399997, 0.078261563764200001, 0.27203219361999997), (-0.042804523999799997, 0.29394584933500001, 0.041995361459999998), (-0.168625231643, -0.035457039764399997, -0.24557754291600001), (-0.19457979239500001, -0.103406961977, 0.20358218145599999), (-0.14621148597600001, 0.23003683216000001, 0.12532061769200001), (-0.19457979239500001, 0.103406961977, -0.20358218145599999), (-0.042804523999799997, -0.29394584933500001, -0.041995361459999998), (0.168625231643, -0.035457039764399997, 0.24557754291600001), (0.19457979239500001, -0.103406961977, -0.20358218145599999), (0.042804523999799997, 0.29394584933500001, -0.041995361459999998), (0.19457979239500001, 0.103406961977, 0.20358218145599999), (0.27203219361999997, -0.099366056939399997, -0.078261563764200001), (-0.099366056939399997, 0.078261563764200001, -0.27203219361999997), (0.042804523999799997, -0.29394584933500001, 0.041995361459999998), (-0.29394584933500001, 0.041995361459999998, 0.042804523999799997), (0.035457039764399997, -0.24557754291600001, 0.168625231643), (-0.23003683216000001, 0.12532061769200001, 0.14621148597600001), (-0.14621148597600001, -0.23003683216000001, -0.12532061769200001), (0.12532061769200001, -0.14621148597600001, 0.23003683216000001), (-0.168625231643, 0.035457039764399997, 0.24557754291600001), (0.078261563764200001, -0.27203219361999997, -0.099366056939399997), (0.27203219361999997, 0.099366056939399997, 0.078261563764200001), (0.041995361459999998, -0.042804523999799997, 0.29394584933500001), (0.168625231643, 0.035457039764399997, -0.24557754291600001), (-0.24557754291600001, -0.168625231643, -0.035457039764399997), (0.14621148597600001, -0.23003683216000001, 0.12532061769200001), (-0.12532061769200001, -0.14621148597600001, -0.23003683216000001), (0.12532061769200001, 0.14621148597600001, -0.23003683216000001), (-0.23003683216000001, -0.12532061769200001, -0.14621148597600001), (-0.24557754291600001, 0.168625231643, 0.035457039764399997), (0.035457039764399997, 0.24557754291600001, -0.168625231643), (0.103406961977, -0.20358218145599999, -0.19457979239500001), (-0.29394584933500001, -0.041995361459999998, -0.042804523999799997), (-0.103406961977, 0.20358218145599999, -0.19457979239500001), (0.29394584933500001, 0.041995361459999998, -0.042804523999799997), (-0.035457039764399997, -0.24557754291600001, -0.168625231643), (0.23003683216000001, 0.12532061769200001, -0.14621148597600001), (-0.27203219361999997, 0.099366056939399997, -0.078261563764200001), (0.099366056939399997, -0.078261563764200001, -0.27203219361999997), (0.20358218145599999, -0.19457979239500001, -0.103406961977), (-0.20358218145599999, -0.19457979239500001, 0.103406961977), (-0.078261563764200001, -0.27203219361999997, 0.099366056939399997), (0.078261563764200001, 0.27203219361999997, 0.099366056939399997), (-0.041995361459999998, -0.042804523999799997, -0.29394584933500001), (0.24557754291600001, -0.168625231643, 0.035457039764399997)]
    facetPointsObjectList= []
    for x,y,z in facetPoints:
        p = Point3D(x,y,z)
        facetPointsObjectList.append(p)
    return facetPointsObjectList

def readInResidueDescritors():
    desciptors= [('NAME','CHARGED','HYD'),('ALA',0,1.8),('ARG',1,-4.5),('ASN',0,-3.5),('ASP',1,-3.5),('CYS',0,2.5),('GLN',0,-3.5),('GLU',1,-3.5),('GLY',0,-0.4),('HIS',0,-3.2),('ILE',0,4.5),('LEU',0,3.8),('LYS',1,-3.9),('MET',0,1.9),('PHE',0,2.8),('PRO',0,-1.6),('SER',0,-0.8),('THR',0,-0.7),('TRP',0,-0.9),('TYR',0,-1.3),('VAL',0,4.2)]
    residueDescriptors= {}
    for name,chgd,hyd in desciptors[1:]:
        residueDescriptors[name]= ResidueDescriptors(chgd,hyd)
    return residueDescriptors

def loadMol(loadMol):
    ifs = oemolistream(loadMol)
    mol = OEGraphMol()
    OEReadMolecule(ifs, mol)
    return mol

def loadLigs(loadMol):
    ifs = oemolistream(loadMol)
    mainMol= OEGraphMol()
    molNum= 0
    for mol in ifs.GetOEGraphMols():
        mol = OEGraphMol(mol)
        for atom in mol.GetAtoms():
            atom.SetType(str(molNum))
        OEAddMols(mainMol,mol)
        molNum += 1
    if molNum >= 3:
        mainMol=trimLigs(mainMol)
    else:
        print "Not enough molecules to trim ligands."
    mainMol.SetTitle('ligs')
    if mainMol.NumAtoms() == 0:
        raise Exception,"No Atoms in ligand"
    return mainMol

def trimLigs(mol):
    cent = OEFloatArray(3)
    ext = OEFloatArray(3)
    OEGetCenterAndExtents(mol,cent,ext)
    grid = OEScalarGrid()
    atomsInCell= {}
    molsInCell= {}
    OEMakeGridFromCenterAndExtents(grid,cent,ext,2.)
    for atom in mol.GetAtoms():
        x,y,z = mol.GetCoords(atom)
        spIdx= grid.SpatialCoordToElement(x,y,z)
        try:
            molsInCell[spIdx].add(atom.GetType())
            atomsInCell[spIdx].append(atom)
        except:
            atomsInCell[spIdx]=[atom]
            molsInCell[spIdx]= set([atom.GetType()])
    for cellIdx,atomIdxes in atomsInCell.items():
        if len(molsInCell[cellIdx]) < 3:
            for atomIdx in atomIdxes:
                mol.DeleteAtom(atomIdx)
    return mol

class Lid:
    def __init__(self,surf,lidSurfPoints,burial,area,centP):
        self.centP= centP
        self.burial= burial
        self.surf= surf
        self.area= area
        self.lidSurfPoints=lidSurfPoints

class Point3D:
    def __init__ (self,x,y,z):
        self.x= float(x)
        self.y= float(y)
        self.z= float(z)

    def distanceSq (self,o):
        if type(o) is tuple:
            x,y,z= o
        else:
            x= o.x
            y= o.y
            z= o.z
        dx= (self.x-x)
        dy= (self.y-y)
        dz= (self.z-z)
        return(dx*dx+dy*dy+dz*dz)

    def angle(self,o):
        pnx,pny,pnz = (self.x,self.y,self.z)
        lnx,lny,lnz = (o.x,o.y,o.z)
        dot= pnx*lnx + pny*lny + pnz*lnz
        dot= min(dot,1.0)
        dot= max(dot,-1.0)
        angle = math.acos(dot)*180./math.pi
        return(angle)

    def angPointToPlane(self,o):
        ax,ay,az = (self.x,self.y,self.z)
        px,py,pz = (o.loc.x,o.loc.y,o.loc.z)
        pnx,pny,pnz = (o.normal.x,o.normal.y,o.normal.z)
        side = (pnx*(ax-px))+(pny*(ay-py))+(pnz*(az-pz))
        return side
    
class ResidueDescriptors:
    def __init__(self,charged,hyd):
        self.charged= float(charged)
        self.hyd= float(hyd)
        self.residueNumbers= set()
        
    def addRes(self,i):
        self.residueNumbers.add(i)
                
    def resCount(self):
        return float(len(self.residueNumbers))
    
class Atom:
    def __init__ (self,x,y,z,atom=None,idx=None,residue=None,hetAtm=None,isHeavy=None, isBackbone=None,resKey=None):
        self.loc= Point3D(x,y,z)
        self.idx = idx
        self.resKey = resKey
        self.atom = atom
        self.x= x
        self.y= y
        self.z= z
        self.isBackbone = isBackbone
        self.isInPocket = False
        self.isHetAtm = hetAtm
        self.loc = Point3D(x,y,z)
        self.residue = residue
        self.isInInitPocket = False
        self.isHeavy = isHeavy

    def angPointToPlane(self,o):
        return(self.loc.angPointToPlane(o))

    def distanceSq (self, o):
        if type(o) is tuple: return(self.loc.distanceSq(o))
        else: return(self.loc.distanceSq(o.loc))

    def x (self):
        return(self.loc.x)

    def y (self):
        return(self.loc.y)

    def z (self):
        return(self.loc.z)

    def march(self,facetPoint,marchDist):
        ax,ay,az= facetPoint.x+self.x,facetPoint.y+self.y,facetPoint.z+self.z
        dist= 0.3
        t= (dist+marchDist)/dist
        cx= self.x+(ax-self.x)*t
        cy= self.y+(ay-self.y)*t
        cz= self.z+(az-self.z)*t
        return (cx,cy,cz)

class SurfacePoint:
    def __init__ (self,x,y,z,nx,ny,nz,idx, inSurf):
        self.loc= Point3D(x,y,z)
        self.normal= Point3D(nx,ny,nz)
        self.neighbors=[]
        self.idx= idx
        self.inInitSet= False
        self.component= -1
        self.initComponent= -1
        self.surf= inSurf
        self.checkedLid= False

    def __str__ (self):
        s= "X= "+str(self.x())+" Y= "+str(self.y())+" Z="+str(self.z())
        return(s)

    def hasNays (self, inSet):
        n= False
        for x in self.neighbors:
            p= self.surf.getPoint(x)
            if p.inInitSet == inSet:
                n= True
        return(n)

    def setInInitSet(self):
        self.inInitSet= True

    def removeInitSet(self):
        self.inInitSet= False

    def distanceSq (self, o):
        return(self.loc.distanceSq(o.loc))

    def angle(self,o):
        return(self.normal.angle(o.normal))

    def setNay (self,n):
        self.neighbors.append(n)

    def x (self):
        return(self.loc.x)

    def y (self):
        return(self.loc.y)

    def z (self):
        return(self.loc.z)

class Surface:
    def __init__ (self,inMol,protAtoms=None,facetPoints=None):
        self.depth= None
        self.facetPoints= facetPoints
        self.sitePointsMol= OEGraphMol()
        self.sitePointsDict= {}
        self.burialLid= OESurface()
        self.lidArea= 0
        self.protAtoms = protAtoms
        self.depthMol = OEGraphMol()
        self.surfAtmDict = {}
        self.residueList = []
        self.oeSurf = OESurface()
        self.oeMol= OEGraphMol(inMol)
        self.surfPoints=[]
        self.grid= OEScalarGrid()
        self.gridHighRes = OEScalarGrid()
        self.gridIdxToSurfaceIdx={}
        self.gridHighResIdxToSurfaceIdx={}
        self.burialDict={}
        OEAssignBondiVdWRadii(self.oeMol)
        OEMakeMolecularSurface(self.oeSurf, self.oeMol,0.5,1.05)
        OESurfaceToMoleculeDistance(self.oeSurf, self.oeMol)
        for i in xrange(self.oeSurf.GetNumVertices()):
            x= self.oeSurf.GetVerticesElement(i * 3)
            y= self.oeSurf.GetVerticesElement(i * 3 + 1)
            z= self.oeSurf.GetVerticesElement(i * 3 + 2)
            nx= self.oeSurf.GetNormalsElement(i * 3)
            ny= self.oeSurf.GetNormalsElement(i * 3 + 1)
            nz= self.oeSurf.GetNormalsElement(i * 3 + 2)
            v= SurfacePoint(x,y,z,nx,ny,nz,i,self)
            self.surfPoints.append(v)
        for i in xrange(self.oeSurf.GetNumTriangles()):
            v1 = self.oeSurf.GetTrianglesElement(i*3)
            v2 = self.oeSurf.GetTrianglesElement(i*3+1)
            v3 = self.oeSurf.GetTrianglesElement(i*3+2)
            self.makeNays(v1,v2)
            self.makeNays(v1,v3)
            self.makeNays(v2,v3)

    def refinePocket(self):
        self.iteratePocket('shrink', 5)
        self.iteratePocket('grow', 15)
        self.iteratePocket('shrink', 10)        
    
    def iteratePocket(self,type,count):
        if type == 'grow':
            inInitSet= False
        else:
            inInitSet= True
        num= 0
        nodesToVisit= []
        while num < count:
            for p in self.surfPoints:
                if p.inInitSet == inInitSet and p.hasNays(not inInitSet):
                    nodesToVisit.append(p)
            while len(nodesToVisit) > 0:
                n= nodesToVisit.pop()
                n.inInitSet= not inInitSet
            num+=1       

    def marchOutCavity(self):
        initNodes=[]
        nodesToVisit= []
        self.iteratePocket('shrink', 1)
        for p in self.surfPoints:
            if p.inInitSet == True and p.hasNays(False):
                initNodes.append(p)
        for m in initNodes:
            visited= []
            nodesToVisit= []
            for nay in m.neighbors:
                nodesToVisit.append(self.getPoint(nay))
            while len(nodesToVisit) > 0:
                n= nodesToVisit.pop()
                visited.append(n.idx)
                if n.angle(m) < 35.:
                    n.inInitSet= True
                    for nay2 in n.neighbors:
                        n2= self.getPoint(nay2)
                        if not n2.idx in visited:
                            nodesToVisit.append(n2)


    def findConnectedGroupsInInitSet(self):
        #Finds connected components in initial set
        self.initComps = {}
        nodesToVisit=[]
        curInitComp= 0
        addedInit= True
        while addedInit:
            curInitComp+=1
            addedInit= False
            for p in self.surfPoints:
                if p.inInitSet == True and p.initComponent == -1 and len(nodesToVisit) == 0:
                    nodesToVisit.append(p)
                    addedInit= True
            while len(nodesToVisit) > 0:
                n= nodesToVisit.pop()
                if n.initComponent == -1:
                    n.initComponent= curInitComp
                    for nay in n.neighbors:
                        sp= self.getPoint(nay)
                        if sp.initComponent == -1 and sp.inInitSet == True:
                            nodesToVisit.append(sp)
        for p  in self.surfPoints:
            if p.initComponent != -1:
                try:
                    self.initComps[p.initComponent].append(p)
                except:
                    self.initComps[p.initComponent]=[]
                    self.initComps[p.initComponent].append(p)
        for x in self.initComps.keys():
            for y in self.initComps.keys():
                if len(self.initComps[x]) > len(self.initComps[y]):
                    for point in self.initComps[y]:
                        point.removeInitSet()
                elif len(self.initComps[x]) < len(self.initComps[y]):
                    for point in self.initComps[x]:
                        point.removeInitSet()
        for p in self.surfPoints:
            p.initComponent = -1

    def makeNays (self,id1,id2):
        self.getPoint(id1).setNay(id2)
        self.getPoint(id2).setNay(id1)

    def getPoint(self,i):
        return(self.surfPoints[i])

    def makeGrid(self,res):
        ctr = OEFloatArray(3)
        ext = OEFloatArray(3)
        OEGetSurfaceCenterAndExtents(self.oeSurf,ctr,ext)
        OEMakeGridFromCenterAndExtents(self.grid,ctr,ext,res)
        for sp in self.surfPoints:
            spIdx = self.grid.SpatialCoordToElement(sp.x(),sp.y(),sp.z())
            try:
                self.gridIdxToSurfaceIdx[spIdx].append(sp.idx)
            except:
                self.gridIdxToSurfaceIdx[spIdx]=[]
                self.gridIdxToSurfaceIdx[spIdx].append(sp.idx)
        # High Resolution Grid for burial determination
        OEMakeGridFromCenterAndExtents(self.gridHighRes,ctr,ext,.5)
        for sp in self.surfPoints:
            spIdx = self.gridHighRes.SpatialCoordToElement(sp.x(),sp.y(),sp.z())
            try:
                self.gridHighResIdxToSurfaceIdx[spIdx].append(sp.idx)
            except:
                self.gridHighResIdxToSurfaceIdx[spIdx]=[]
                self.gridHighResIdxToSurfaceIdx[spIdx].append(sp.idx)

    def makeAtomSet(self,protAtoms):
        resDict= {}
        resCount= {}
        for sp in self.surfPoints:
            if sp.inInitSet:
                atom = self.oeMol.GetAtom(OEHasAtomIdx(self.oeSurf.GetAtomsElement(sp.idx)))
                residue = OEAtomGetResidue(atom)
                resNum = residue.GetResidueNumber()
                try:
                    resDict[resNum].append(atom)
                except:
                    resDict[resNum]= [atom]
                resCount[resNum]= resCount.get(resNum,0) + 1
        for k,v in resCount.items():
            if v > 2:
                for atom in resDict[k]:
                    protAtoms[atom.GetIdx()].isInInitPocket= True
                    protAtoms[atom.GetIdx()].isInPocket= True


    def makeColorSurface(self,isPocket=False):
        self.nSurf= OESurface(self.oeSurf)
        for sp in self.surfPoints:
            if sp.inInitSet:
                self.oeSurf.SetVertexCliqueElement(sp.idx, 12)
                self.nSurf.SetColorElement(sp.idx, 1.,0.,0.,.6)
                self.nSurf.SetVertexCliqueElement(sp.idx, 12)
            else:
                self.oeSurf.SetColorElement(sp.idx, 1.,0.,0.,.6) 
                self.nSurf.SetColorElement(sp.idx, 1.,0.,0.,.6) 
        OESurfaceCropToClique(self.nSurf,12)


    def fillActiveSite(self):
        for p in self.surfPoints:
            if p.inInitSet == True:
                self.oeSurf.SetVertexCliqueElement(p.idx, 5)
            else:
                self.oeSurf.SetVertexCliqueElement(p.idx, 1)                 
        self.cropSurf = OESurface(self.oeSurf)
        OESurfaceCropToClique(self.cropSurf,5)
        self.cropGrid = OEScalarGrid()
        ctr = OEFloatArray(3)
        ext = OEFloatArray(3)
        OEGetSurfaceCenterAndExtents(self.cropSurf,ctr,ext)
        OEMakeGridFromCenterAndExtents(self.cropGrid,ctr,ext,.7) 
        OEMakeGridFromSurface(self.cropGrid,self.cropSurf)
        for ii in xrange(0,self.cropGrid.GetXDim()):
            for ij in xrange(0,self.cropGrid.GetYDim()):
                for ik in xrange(0,self.cropGrid.GetZDim()):
                    dist = self.cropGrid.GetValue(ii,ij,ik)
                    if dist < 64.:
                        (x,y,z)= self.cropGrid.GridIdxToSpatialCoord(ii,ij,ik)
                        (bi,bj,bk)= self.grid.SpatialCoordToGridIdx(x,y,z)
                        minD = (65.,None)
                        minx=  max(0,bi-8)
                        miny=  max(0,bj-8)
                        minz=  max(0,bk-8)
                        maxx=  min(self.grid.GetXDim(),bi+8)
                        maxy=  min(self.grid.GetYDim(),bj+8)
                        maxz=  min(self.grid.GetZDim(),bk+8)
                        for k in xrange(minx, maxx):
                            for j in xrange(miny, maxy):
                                for i in xrange(minz, maxz):
                                    idx= self.grid.GridIdxToElement(k,j,i)  
                                    if self.gridIdxToSurfaceIdx.get(idx):                           
                                        for spIdx in self.gridIdxToSurfaceIdx[idx]:
                                            p= self.getPoint(spIdx)
                                            if p.inInitSet:
                                                dist2 = ((p.loc.x-x)**2+(p.loc.y-y)**2+(p.loc.z-z)**2)
                                                if dist2 < minD[0]:
                                                    minD = (dist2,p)

                        if minD[0] < 64. and minD[0] > 1.: #0.1225:
                            minPoint = minD[1]
                            point= Point3D(x,y,z)
                            side= point.angPointToPlane(minPoint)
                            if side >= 0.:
                                floatArray = OEFloatArray([x,y,z])
                                newIdx = self.sitePointsMol.NewAtom(1).GetIdx()
                                self.sitePointsMol.GetAtom(OEHasAtomIdx(newIdx))
                                self.sitePointsDict[newIdx]= Atom(x,y,z,self.sitePointsMol.GetAtom(OEHasAtomIdx(newIdx)),newIdx)
                                self.sitePointsMol.SetCoords(self.sitePointsMol.GetAtom(OEHasAtomIdx(newIdx)),floatArray)
                                self.sitePointsMol.GetAtom(OEHasAtomIdx(newIdx)).SetName('sphere')


    def pointBurial(self,protSurf,p,ignoreInPocket=False):
        facetCount= 0
        for facetPoint in self.facetPoints:
            marchDist=0
            noCount= True
            while marchDist <= 50. and noCount:
                mX,mY,mz= p.march(facetPoint,marchDist)
                try:
                    gElement= protSurf.gridHighRes.SpatialCoordToElement(mX,mY,mz)
                except:
                    gElement= False
                if gElement:             
                    if protSurf.gridHighResIdxToSurfaceIdx.get(gElement):
                        for spIdx in protSurf.gridHighResIdxToSurfaceIdx[gElement]:
                            pSurfPoint= protSurf.getPoint(spIdx)
                            if pSurfPoint.inInitSet or ignoreInPocket:                               
                                if noCount:
                                    facetCount +=1
                                    noCount= False
                marchDist +=.25
        return facetCount

    def checkSitePointBurial(self):
        for sitePoint in self.sitePointsDict.values():
            facetCount= self.pointBurial(self,sitePoint)
            if facetCount < 40:
                self.sitePointsMol.DeleteAtom(self.sitePointsMol.GetAtom(OEHasAtomIdx(sitePoint.idx)))
            else:
                self.burialDict[sitePoint.idx]=facetCount
        multipleCliqs= True
        while multipleCliqs:
            multipleCliqs= self.removeDisconnectedSitePoints()

    def checkPocketSelection(self):
        self.gridIdxToSitePoints= {}
        for a in self.sitePointsMol.GetAtoms():
            x,y,z = self.sitePointsMol.GetCoords(a)
            spIdx= self.grid.SpatialCoordToElement(x,y,z)
            try:
                self.gridIdxToSitePoints[spIdx].append((x,y,z))
            except:
                self.gridIdxToSitePoints[spIdx]=[(x,y,z)]
        for sp in self.surfPoints:
            if sp.inInitSet:
                minD = 1000.
                (bi,bj,bk)= self.grid.SpatialCoordToGridIdx(sp.x(),sp.y(),sp.z())
                minx=  max(0,bi-4)
                miny=  max(0,bj-4)
                minz=  max(0,bk-4)
                maxx=  min(self.grid.GetXDim(),bi+4)
                maxy=  min(self.grid.GetYDim(),bj+4)
                maxz=  min(self.grid.GetZDim(),bk+4)
                while minD > 9:
                    for k in xrange(minx, maxx):
                        for j in xrange(miny, maxy):
                            for i in xrange(minz, maxz):
                                idx= self.grid.GridIdxToElement(k,j,i)
                                try:
                                    siteCoordsList= self.gridIdxToSitePoints[idx]
                                except:
                                    siteCoordsList= None
                                if siteCoordsList:
                                    for x,y,z in siteCoordsList:
                                        dist2 = ((sp.x()-x)**2+(sp.y()-y)**2+(sp.z()-z)**2)
                                        minD = min(minD,dist2)
                    break
                if minD > 9:
                    sp.inInitSet= False
        contGoing= True

        while contGoing:
            contGoing= False
            hasPoints= False
            for sp in self.surfPoints:
                if sp.inInitSet and sp.hasNays(False):
                    hasPoints= True
                    minD = 1000.
                    (bi,bj,bk)= self.grid.SpatialCoordToGridIdx(sp.x(),sp.y(),sp.z())
                    minx=  max(0,bi-2)
                    miny=  max(0,bj-2)
                    minz=  max(0,bk-2)
                    maxx=  min(self.grid.GetXDim(),bi+2)
                    maxy=  min(self.grid.GetYDim(),bj+2)
                    maxz=  min(self.grid.GetZDim(),bk+2)
                    while minD > 3.61:
                        for k in xrange(minx, maxx):
                            for j in xrange(miny, maxy):
                                for i in xrange(minz, maxz):
                                    idx= self.grid.GridIdxToElement(k,j,i)
                                    try:
                                        siteCoordsList= self.gridIdxToSitePoints[idx]
                                    except:
                                        siteCoordsList= None
                                    if siteCoordsList:
                                        for x,y,z in siteCoordsList:
                                            dist2 = ((sp.x()-x)**2+(sp.y()-y)**2+(sp.z()-z)**2)
                                            minD = min(minD,dist2)
                        break
                    if minD > 3.61:
                        sp.inInitSet= False
                        contGoing= True
            if not hasPoints:
                break

    def checkEnclosedPocket(self,protSurf):
        wholePocket = True
        for sp in protSurf.surfPoints:
            if sp.inInitSet and sp.hasNays(False):
                wholePocket= False
        return wholePocket

    def removeDisconnectedSitePoints(self):
        deletedAnAtom= False
        self.tightSurface= OESurface()
        OEAssignBondiVdWRadii(self.sitePointsMol)
        OEMakeMolecularSurface(self.tightSurface, self.sitePointsMol,0.05,0.1)
        nclqs = OEMakeConnectedSurfaceCliques(self.tightSurface)
        areas = []
        for i in range(1, nclqs+1):
            areas.append((OESurfaceCliqueArea(self.tightSurface, i), i))
        areas.sort()
        badCliques= [c for a,c in areas[:-1]]
        for i in xrange(self.tightSurface.GetNumVertices()):
            clique= self.tightSurface.GetVertexCliqueElement(i)
            if clique in badCliques:
                atom = self.sitePointsMol.GetAtom(OEHasAtomIdx(self.tightSurface.GetAtomsElement(i)))
                try:
                    del self.burialDict[atom.GetIdx()]
                    deletedAnAtom= True
                    self.sitePointsMol.DeleteAtom(atom)
                except:
                    pass
        if len(areas) > 1 and deletedAnAtom:
            multipleCliqs= True
        else:
            multipleCliqs= False
        return multipleCliqs

    def noLid(self,protSurf):
        for sp in protSurf.surfPoints:
            if sp.inInitSet:
                protSurf.oeSurf.SetVertexCliqueElement(sp.idx, 12)
        self.sphereSurf= OESurface(protSurf.oeSurf)
        OESurfaceCropToClique(self.sphereSurf,12)
        if self.checkEnclosedPocket(protSurf):
            print 'enclosed'
            maxDist = (0,None,None)
            spCopy= protSurf.surfPoints[:]
            while len(spCopy) > 1:
                sp= spCopy.pop()
                if sp.inInitSet:
                    for sp2 in spCopy:
                        if sp2.inInitSet:
                            dist = (sp.distanceSq(sp2),sp,sp2)
                            maxDist = max(maxDist,dist)
                        #   print dist, maxDist
        else:
            print 'not enclosed'
            maxDist = (0,None,None)
            for sp in protSurf.surfPoints:
                if sp.inInitSet and sp.hasNays(False):
                    for sp2 in protSurf.surfPoints:
                        if sp2.inInitSet:
                            dist = (sp.distanceSq(sp2),sp,sp2)
                            maxDist = max(maxDist,dist)
        lidFloatArray = OEFloatArray([maxDist[1].loc.x,maxDist[1].loc.y,maxDist[1].loc.z])
        pocketFloatArray = OEFloatArray([maxDist[2].loc.x,maxDist[2].loc.y,maxDist[2].loc.z])
        self.makeDepthMol(lidFloatArray,pocketFloatArray)
    
    def makeDepthMol(self,p1,p2):
        self.depthMol.SetCoords(self.depthMol.NewAtom(75),p1)
        self.depthMol.SetCoords(self.depthMol.NewAtom(75),p2)
        self.depth= OEGetDistance(self.depthMol, self.depthMol.GetAtom(OEHasAtomIdx(0)), self.depthMol.GetAtom(OEHasAtomIdx(1)))
               
    def findLid(self):
        self.sphereSurf= OESurface(self.oeSurf)
        self.mainLidSurf= OESurface(self.oeSurf)
        initSetList= []
        for sp in self.surfPoints:
            atom = self.oeMol.GetAtom(OEHasAtomIdx(self.sphereSurf.GetAtomsElement(sp.idx)))
            if atom.GetName() == 'sphere':
                sp.inInitSet= True
                initSetList.append(sp)
        return initSetList

    def findLidCenter(self,potentialLid):
        lidSurfPoints = []
        for i in xrange(potentialLid.GetNumVertices()):
            x= potentialLid.GetVerticesElement(i * 3)
            y= potentialLid.GetVerticesElement(i * 3 + 1)
            z= potentialLid.GetVerticesElement(i * 3 + 2)
            nx= potentialLid.GetNormalsElement(i * 3)
            ny= potentialLid.GetNormalsElement(i * 3 + 1)
            nz= potentialLid.GetNormalsElement(i * 3 + 2)
            v= SurfacePoint(x,y,z,nx,ny,nz,i,self)
            lidSurfPoints.append(v)
        xAveLid= numpy.mean([a.x() for a in lidSurfPoints])
        yAveLid= numpy.mean([a.y() for a in lidSurfPoints])
        zAveLid= numpy.mean([a.z() for a in lidSurfPoints])
        centP= Atom(xAveLid,yAveLid,zAveLid)
        lid_center= (1000,None)
        for point in lidSurfPoints:
            dist= centP.distanceSq(point)
            lid_center= min(lid_center,(dist,point.loc))
        center_point= (10000,None)
        for atom in self.oeMol.GetAtoms():
            if atom.GetName() == 'sphere':
                x,y,z= self.oeMol.GetCoords(atom)
                dist= lid_center[1].distanceSq((x,y,z))
                center_point= min(center_point,(dist,(x,y,z)))
        return Atom(center_point[1][0],center_point[1][1],center_point[1][2]),lidSurfPoints
                
    def findMainLid(self):
        self.lids.sort(key=operator.attrgetter('area'), reverse=True)
        self.mainLidSurfPoints= self.lids[0].lidSurfPoints
        self.mainLidSurf= self.lids[0].surf
        self.mainLidCenter= self.lids[0].centP
        if len(self.lids) > 1 and self.lids[0].burial > 55:
            for smallerLid in self.lids[1:]:
                if smallerLid.area/self.lids[0].area > 0.6:
                    if smallerLid.burial < self.lids[0].burial:
                        self.mainLidSurfPoints= smallerLid.lidSurfPoints
                        self.mainLidSurf= smallerLid.surf
                        self.mainLidCenter= smallerLid.centP
                        break                        
                 
    def makeSphereSurf(self,protSurf=None):
        lidPoints= self.findLid()
        if not lidPoints:
            self.noLid(protSurf)
        else:
            for sl in self.surfPoints:
                if sl.inInitSet:
                    self.sphereSurf.SetVertexCliqueElement(sl.idx, 15)
            OESurfaceCropToClique(self.sphereSurf, 15)
            self.burialLid= OESurface(self.sphereSurf)
            nclqs = OEMakeConnectedSurfaceCliques(self.sphereSurf)
            self.lids = []
            for clq in range(1, nclqs+1):
                potentialLid= OESurface()
                OEMakeCliqueSurface(potentialLid,self.sphereSurf,clq)
                area= OESurfaceCliqueArea(self.sphereSurf, clq)
                self.lidArea += area
                centP,lidSurfPoints= self.findLidCenter(potentialLid)
                facetCount= self.pointBurial(protSurf,centP, True)
                
                self.lids.append(Lid(potentialLid,lidSurfPoints,facetCount,area,centP))
            self.findMainLid()
            maxDist = (0,None)
            for pocketPoint in protSurf.surfPoints:
                minDist = (10000,None)
                if pocketPoint.inInitSet:
                    for lidPoint in self.mainLidSurfPoints:
                        dist = (pocketPoint.distanceSq(lidPoint),pocketPoint)
                        minDist = min(dist,minDist)
                    maxDist = max(maxDist,minDist)
            lidFloatArray = OEFloatArray([self.mainLidCenter.loc.x,self.mainLidCenter.loc.y,self.mainLidCenter.loc.z])
            pocketFloatArray = OEFloatArray([maxDist[1].loc.x,maxDist[1].loc.y,maxDist[1].loc.z])         
            self.makeDepthMol(lidFloatArray, pocketFloatArray)

    def iterateGrid(self,surf,sp,maxDist):
        matchedSurfacePoints= []
        try:
            (ix,iy,iz)= surf.grid.SpatialCoordToGridIdx(sp.x(),sp.y(),sp.z())
            inGrid= True
        except:
            inGrid= False
        if inGrid:
            minx=  max(0,ix-maxDist)
            miny=  max(0,iy-maxDist)
            minz=  max(0,iz-maxDist)
            maxx=  min(surf.grid.GetXDim(),ix+maxDist)
            maxy=  min(surf.grid.GetYDim(),iy+maxDist)
            maxz=  min(surf.grid.GetZDim(),iz+maxDist)
            for k in xrange(minx, maxx):
                for j in xrange(miny, maxy):
                    for i in xrange(minz, maxz):
                            idx= surf.grid.GridIdxToElement(k,j,i)
                            spIdxes= surf.gridIdxToSurfaceIdx.get(idx)
                            if spIdxes:
                                for spIdx in spIdxes:                              
                                    matchedSurfacePoints.append(surf.getPoint(spIdx))
        return matchedSurfacePoints
    
    
    def intialProteinPointsCloseToLigand(self,ligSurface,dTol,angTol):
        dSqTol= float(dTol**2)
        for sp in ligSurface.surfPoints:
            points= self.iterateGrid(self,sp,dTol+1)
            for pSurfPoint in points:
                dSq= sp.distanceSq(pSurfPoint)
                if dSq < dSqTol:
                    ang= sp.angle(pSurfPoint)
                    if ang > angTol:
                        pSurfPoint.setInInitSet()

class AnalyzeProtein:
    def __init__(self, protein,ligand,facetPoints,dTol,angTol,protname,backboneAtoms,residueDescriptors):
        self.residueDescriptors= residueDescriptors
        self.protname= protname
        self.dTol= dTol
        self.angTol= angTol
        self.backboneAtoms= backboneAtoms
        self.facetPoints = facetPoints
        self.protein = protein
        self.ligand = ligand
        self.surfProt= OEGraphMol(self.protein)
        self.protAtoms = self.loadProtAtoms(self.protein)
        self.surfProtAtoms = self.loadProtAtoms(self.surfProt,stripHet=True)
        self.ligSurface =  Surface(self.ligand)
        self.protSurface = Surface(self.surfProt,self.surfProtAtoms,facetPoints=self.facetPoints)
        
    def finalizePocket(self):
        self.pocketMol= OEGraphMol(self.surfProt)
        self.pocketMolWithSitePoints= OEGraphMol(self.pocketMol)
        OEAddMols(self.pocketMolWithSitePoints,self.protSurface.sitePointsMol)
        self.pocketMolWithSitePointsSurf= OESurface()
        OEMakeMolecularSurface(self.pocketMolWithSitePointsSurf,self.pocketMolWithSitePoints)
        self.addHetAtms()
        self.pCharged,self.aHyd= self.calculateOneDimensionalDescriptiors(self.pocketMol)
            
    def makePocket(self):
        self.protSurface.makeGrid(1.)
        self.protSurface.intialProteinPointsCloseToLigand(self.ligSurface,self.dTol,self.angTol)     
        self.protSurface.marchOutCavity()
        self.protSurface.refinePocket()
        self.protSurface.findConnectedGroupsInInitSet()

    def analyzePocket(self):
        self.protSurface.fillActiveSite()
        self.protSurface.checkSitePointBurial()
        self.protSurface.checkPocketSelection()
        self.protSurface.makeColorSurface()
        self.volumePocket= Surface(self.protSurface.sitePointsMol,facetPoints=self.facetPoints)
        self.volume= OESurfaceVolume(self.volumePocket.oeSurf)
        self.hybridMol= OEGraphMol(self.protSurface.sitePointsMol)
        OEAddMols(self.hybridMol,self.protSurface.oeMol)
        self.hybridSurf= Surface(self.hybridMol,facetPoints=self.facetPoints)
        self.hybridSurf.makeSphereSurf(protSurf=self.protSurface)
        self.makePocketMol(self.protSurface,self.surfProtAtoms)
        self.area= OESurfaceCliqueArea(self.protSurface.nSurf,12)
        self.degreeOfBurial= (self.area/(self.area+self.hybridSurf.lidArea))*2.-1.
    
    def writeMol(self):
        pofs = oemolostream(self.protname.split('_')[0]+"_pocket.pdb")
        oebOfs = oemolostream(str(self.protname)+"_all.oeb")
        OESetSurfaceColor(self.hybridSurf.mainLidSurf,0.,1.,1.,.6)
        self.pocketMol.SetData("Surface",self.hybridSurf.burialLid)
        self.hybridSurf.depthMol.SetData("Surface",self.hybridSurf.mainLidSurf)
        OESetSurfaceColor(self.volumePocket.oeSurf,1.,0.,1.,.8)
        self.volumePocket.oeMol.SetData("Surface", self.volumePocket.oeSurf)
        self.volumePocket.oeMol.SetTitle("Volume")
        OEWriteMolecule(oebOfs,self.volumePocket.oeMol)
        self.hybridSurf.depthMol.SetTitle('Pocket_Depth')
        self.ligand.SetTitle("Ligand")
        self.protein.SetTitle("Protein")
        self.protein.SetData("Surface", self.protSurface.nSurf)
        self.pocketMol.SetTitle('Pocket')
        OEWriteMolecule(oebOfs,self.hybridSurf.depthMol)
        OEWriteMolecule(oebOfs,self.ligand)
        OEWriteMolecule(oebOfs,self.protein)
        OEWriteMolecule(oebOfs,self.pocketMol)
        OEWriteMolecule(pofs,self.pocketMol)
        name = protname.split('_')[0]
        print >>open(name+"_pocket.txt",'w'),"Protein: %s \nVolume: %.2f \nDepth: %.2f \nDegree of Burial: %.2f \nPercent Charged: %.1f \nHydrophobicity: %.2f"%(name,self.volume,self.hybridSurf.depth, self.degreeOfBurial,self.pCharged,self.aHyd)

    def loadLigAtoms(self,lig):
        ligDict = {}
        for atom in lig.GetAtoms():
            cx,cy,cz = lig.GetCoords(atom)
            idx = atom.GetIdx()
            ligDict[atom.GetIdx()] = Atom(cx,cy,cz,atom,idx)
        for atom in lig.GetAtoms(OEIsHeavy()):
            ligDict[atom.GetIdx()].isHeavy = True
        return ligDict

    def loadProtAtoms(self,mol,stripHet=False):
        atomDict = {}
        for atom in mol.GetAtoms():
            atomStripped= False
            hetAtm= False
            isHeavy= False
            isBackbone = False
            cx,cy,cz = mol.GetCoords(atom)
            idx = atom.GetIdx()
            residue = OEAtomGetResidue(atom)
            resNum = residue.GetResidueNumber()
            resChain = residue.GetChainID()
            resKey = str(resNum)+"_"+str(resChain)
            if atom.GetAtomicNum() > 1:
                isHeavy= True
            if residue.IsHetAtom():
                if stripHet:
                    mol.DeleteAtom(atom)
                    atomStripped= True
                else:
                    hetAtm = True
            if not atomStripped:
                if atom.GetName() in self.backboneAtoms:
                    isBackbone = True
                atomDict[atom.GetIdx()] = Atom(cx,cy,cz,atom,idx,residue,hetAtm,isHeavy,isBackbone,resKey)
        return atomDict

    def makePocketMol(self,surf,atomDict):
        backBoneList = []
        allAtomList = []
        surf.makeAtomSet(atomDict)
        for atom in atomDict.values():
            if atom.isInInitPocket:
                if not atom.isBackbone:
                    allAtomList.append(atom.resKey)
                else:
                    backBoneList.append(atom.resKey)
        for pAtom in atomDict.values():
            if pAtom.resKey in allAtomList:
                pAtom.isInPocket = True
            elif pAtom.resKey in backBoneList:
                if pAtom.atom.GetName() in self.backboneAtoms:
                    pAtom.isInPocket = True
                    pAtom.residue.SetName('GLY')
        for rmAtom in atomDict.values():
            if not rmAtom.isInPocket:
                self.surfProt.DeleteAtom(rmAtom.atom)

    def addHetAtms(self):
        hetRes = set()
        bitGrid= OEScalarGrid()
        ctr= OEFloatArray(3)
        ext= OEFloatArray(3)
        OEGetSurfaceCenterAndExtents(self.pocketMolWithSitePointsSurf,ctr,ext)
        OEMakeGridFromCenterAndExtents(bitGrid,ctr,ext,.5)
        OEMakeBitGridFromSurface(bitGrid,self.pocketMolWithSitePointsSurf)
        for atom in self.protAtoms.values():
            if atom.isHeavy:
                if atom.isHetAtm:
                    if bitGrid.IsInGrid(atom.x,atom.y,atom.z):
                        (ix,iy,iz)= bitGrid.SpatialCoordToGridIdx(atom.x,atom.y,atom.z)
                        if bitGrid.GetValue(ix,iy,iz) == 1.0:
                            hetRes.add(atom.residue)
        for res in hetRes:
            for hatom in OEGetResidueAtoms(self.protein,res):
                self.pocketMol.NewAtom(hatom)

    def calculateOneDimensionalDescriptiors(self,mol):
        hetAtms= set()
        for atm in mol.GetAtoms():
            res = OEAtomGetResidue(atm)
            resName = res.GetName()
            chain= res.GetChainID()
            resNum= res.GetResidueNumber()
            if res.IsHetAtom() and not atm.GetAtomicNum() in [1,8]:
                hetAtms.add((chain,resNum))
            elif not res.IsHetAtom():
                try:self.residueDescriptors[resName].addRes((chain,resNum))
                except: print 'Unknown Residue: ', resName
        nCharged= float(len(hetAtms))
        totHyd= numRes= 0.
        for resInfo in self.residueDescriptors.values():
            nCharged += (resInfo.charged * resInfo.resCount())
            totHyd += (resInfo.hyd * resInfo.resCount())
            numRes += resInfo.resCount()
        
        return (nCharged/(numRes+float(len(hetAtms))))*100.,(totHyd/numRes)
            
            
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print >> sys.stderr,"usage: %s docked_ligands.sdf protein.pdb" % (sys.argv[0])
        exit(0)
    backboneAtoms= [" N  "," CA "," C  "," O  "," H  "," HA "," N1 "," C1 "," C2 "," O1 "," H1 "," H2 "]
    dTol= 5
    angTol = 170.
    sw = OEStopwatch()
    sw.Start()
    descriptors_list= readInResidueDescritors()
    facetPoints= readInFacetPoints()
    protname = sys.argv[2].split(".")[0].split("/")[-1]
    prot = loadMol(sys.argv[2])
    lig = loadLigs(sys.argv[1])
    prot.SetTitle(sys.argv[2])
    ap= AnalyzeProtein(prot,lig,facetPoints,dTol,angTol,protname,backboneAtoms,descriptors_list)
    ap.makePocket()
    ap.analyzePocket()
    ap.finalizePocket()
    ap.writeMol()
    print "Protein: %s: %.0f Seconds " % (protname,sw.Elapsed())
