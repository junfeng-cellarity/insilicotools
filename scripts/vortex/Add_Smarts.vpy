#Imports of libraries
import java.lang as lang
import java.lang.String as String
import java.net.URL as URL
import java.net.URLConnection as URLConnection
import java.io.InputStreamReader as InputStreamReader
import java.io.BufferedReader as BufferedReader
import java.io.DataOutputStream as DataOutputStream
import java.util.Properties as Properties
#import javax.xml
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.IntArray as IntArray
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.plots.VortexScatterPlot as VortexScatterPlot

__doc__ = 'Add SMARTS Pattern match count as a calculated column'
__author__ = 'Jun Feng'
__version__ = '1.0'

import sys
import time
import java.net.URL
import math
import jarray
import os

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

from java.lang import Runnable, Double
from java.util import Random, BitSet
from java.util.concurrent import ExecutionException

from javax.swing.border import TitledBorder
from javax.swing import BorderFactory, ProgressMonitor, JTextArea, JPanel, JLabel, JCheckBox
from javax.swing import JButton, JFrame, JTextField,JComboBox
from java.awt import GridBagLayout as awtGridBagLayout, Dimension, BorderLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets
import javax.swing.SwingWorker as SwingWorker

import openeye.oechem.OEFormat as OEFormat
import openeye.oechem.OEGraphMol as OEGraphMol
import openeye.oechem.OESubSearch as OESubSearch
import openeye.oechem.OEPartPredAtom as OEPartPredAtom
import openeye.oechem.oechem.OESubsetMol as OESubsetMol
import openeye.oechem.oechem.OEAssignImplicitHydrogens as OEAssignImplicitHydrogens
import openeye.oechem.oechem.OEWriteMolecule as OEWriteMolecule
import openeye.oechem.oechem.OEAddMols as OEAddMols
import openeye.oechem.oechem.OEDetermineComponents as OEDetermineComponents
import openeye.oechem.oechem.OEThrow as OEThrow
import openeye.oechem.oechem.OECreateAbsSmiString as OECreateAbsSmiString
import openeye.oechem.oemolistream as oemolistream
import openeye.oechem.oechem.OEReadMolecule as OEReadMolecule
import openeye.oechem.oechem.OECreateCanSmiString as OECreateCanSmiString
import openeye.oechem.oemolistream as oemolistream


class Smarts:
    def __init__(self):
        self.ringpats = self.initRingPatterns(3,10)

    def initRingPatterns(self,minsize,maxsize):
        if not minsize<=maxsize:
            OEThrow.Fatal('Must have ring minsize<=maxsize.')
        ringsmarts=['' for i in range(minsize)]
        ringsmarts.extend(['*1'+((i-1)*'~*')+'1' for i in range(minsize,maxsize+1)])
        ringpats=[None for i in range(maxsize+1)]
        for size in range(minsize,maxsize+1):
            smarts=ringsmarts[size]
            ringpats[size]=OESubSearch()
            if not ringpats[size].Init(smarts):
                OEThrow.Fatal('Bad SMARTS: %s'%smarts)
        return ringpats

    def rings(self,mol,ringpats,sasr=False,splitrings=False,verbose=False):
        maxsize=len(ringpats)-1
        for i in range(3,maxsize+1):
            if ringpats[i]!=None: break
        minsize=i
        groupmol=OEGraphMol()
        out_count_this=0
        ring_count_this=0
        ringmol=OEGraphMol()
        ir=0
        ringatoms={}
        ringbonds={}
        rsizes=[]

        for size in range(minsize,maxsize+1):
            pat=ringpats[size]
            for matchbase in pat.Match(mol,1):
                if sasr: newring=0
                else: newring=1
                if verbose: str=''
                aids=[]
                bids=[]
                for matchpair in matchbase.GetAtoms():
                    idx = matchpair.target.GetIdx()
                    aids.append(idx)
                    if sasr:
                        if not ringatoms.has_key(idx): ringatoms[idx]=size
                        if ringatoms[idx]==size: newring=1
                    if verbose>1: str+="%s "%matchpair.target.GetName()

                for matchpair in matchbase.GetBonds():
                    idx = matchpair.target.GetIdx()
                    bids.append(idx)
                    if sasr:
                        if not ringbonds.has_key(idx): ringbonds[idx]=size
                        if ringbonds[idx]==size: newring=1

                if newring:
                    ir+=1
                    rlist=[]
                    ringmol.Clear()
                    for i in range(mol.NumAtoms()): rlist.append(0)
                    for j in aids: rlist[j]=1
                    pred = OEPartPredAtom(rlist)
                    pred.SelectPart(1)
                    OESubsetMol(ringmol,mol,pred)
                    OEAssignImplicitHydrogens(ringmol)
                    if splitrings:
                        ringmol.SetTitle(mol.GetTitle()+" (ring)")
                        OEWriteMolecule(oms,ringmol)
                        out_count_this+=1
                    else:
                        OEAddMols(groupmol,ringmol)
                        ring_count_this+=1
                        rsizes.append(size)
                    if verbose:
                        OEThrow.Info("\t%d. ring(%d): %s %s" % (ir,size,str,OECreateCanSmiString(ringmol)))

        return out_count_this,ring_count_this,groupmol,rsizes

    def splitMolecule(self,mol):
        partlist= jarray.zeros(mol.GetMaxAtomIdx(),'i')
        print partlist
        count=OEDetermineComponents(mol,partlist)
        print partlist
        pred = OEPartPredAtom(partlist)
        mol_parts = []
        for i in xrange(1,count+1):
            pred.SelectPart(i)
            partmol = OEGraphMol()
            OESubsetMol(partmol, mol, pred)
            mol_parts.append(partmol)
        return(mol_parts)

    def countAromaticRings(self,mol):
        outCount,ringCount,groupMol,rsizes = self.rings(mol,self.ringpats,sasr=True)
        ringList = self.splitMolecule(groupMol)
        aromaticRingCount = 0
        for ring in ringList:
            ringIsAromatic = True
            for atm in ring.GetAtoms():
                if atm.IsAromatic() == False:
                    ringIsAromatic = False
                    break
            if ringIsAromatic:
                aromaticRingCount += 1
        return aromaticRingCount


class SMARTSFactory:
    def __init__(self):
        self.subsearch = OESubSearch()
        return


    def show_interface(self):
        mypanel = swing.JPanel(BorderLayout())
        mypanel.setBorder(TitledBorder("Input SMARTS"))

        self.textarea = JTextArea()
        self.textarea.setPreferredSize(Dimension(200,200))
        self.textarea.setLineWrap(True)

        mypanel.add(self.textarea,BorderLayout.CENTER)

        self.mycheckbox = swing.JCheckBox("Selected Only")
        structure_cols = []
        for column in vtable.getColumnsArray():
            colName = String(column.getName())
            if colName.contains("SMILES") or colName.contains("smiles") or colName.equalsIgnoreCase("smiles") or colName.equalsIgnoreCase("Structure") or colName.equalsIgnoreCase("MOLFILE"):
                structure_cols.append(colName)
        self.comboBox = JComboBox(structure_cols)

        btnPanel = JPanel()
        self.cb = JCheckBox("Use Count:")
        btnPanel.add(self.cb)
        btnPanel.add(self.comboBox)
        mypanel.add(btnPanel,BorderLayout.SOUTH)

        ret = vortex.showInDialog(mypanel, "Add SMARTS:")

        if ret == vortex.OK:
            text = self.textarea.getText()
            smarts = text.strip()
            return len(smarts)>0
        else:
            return False


    def calculateAll(self):
        smarts = self.textarea.getText().strip()
        subsearch = self.subsearch
        structureCol = self.comboBox.getSelectedItem()
        useCount = self.cb.isSelected()
        class DescriptorTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)

            def doInBackground(self):
                useSelected = False
                nrow = vtable.getRealRowCount()
                rows = range(0,nrow)
                if useSelected:
                    nrow = vtable.getSelectedRowCount()
                    rows = []
                    tmp = vtable.getSelectedRows()
                    for i in range(0,tmp.size()):
                        rows.append(tmp.get(i))

#                items = self.myduallist.getItems()
                if not subsearch.Init(smarts):
                    vortex.alert("Invalid SMARTS pattern.")
                    return

                first = True
                for row in rows:
                    if self.progressMonitor.isCanceled():
                        break

                    self.super__publish([(100*row/len(rows),"process compound %d"%(row))])
                    oemol = OEGraphMol()
                    if String(structureCol).equalsIgnoreCase("Structure"):
                        molfile = vtable.getMolFileManager().getMolFileAtRow(row)
                        molin = oemolistream()
                        molin.SetFormat(OEFormat.MDL)
                        molin.openstring(molfile)
                        OEReadMolecule(molin,oemol)

                    else:
                        structure_column = vtable.findColumnWithName(structureCol,0)
                        smiles = structure_column.getValueAsString(row)
                        molin = oemolistream()
                        molin.SetFormat(OEFormat.SMI)
                        molin.openstring(smiles)
                        OEReadMolecule(molin,oemol)

                    if first:
                        new_column = vtable.findColumnWithName(smarts,1,Column.COLUMN_INT)

                    matchCount = 0
                    for match in subsearch.Match(oemol, 1):
                        matchCount += 1
                        if not useCount:
                            break


                    new_column.setInt(row,matchCount)

                    if first:
                        first = False


            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                self.progressMonitor.close()
                try:
                    self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()
                except ExecutionException, e:
                    e.printStackTrace()
                    vortex.alert(e.getMessage())

        task = DescriptorTask()
        task.execute()

if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        smartsFactory = SMARTSFactory()
        if smartsFactory.show_interface():
            smartsFactory.calculateAll()
