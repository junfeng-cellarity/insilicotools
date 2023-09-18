import sys
import math
import traceback
import os
import javax.swing.ProgressMonitor as ProgressMonitor
import javax.swing.SwingWorker as SwingWorker
import javax.swing.JComboBox as JComboBox
import javax.swing.JTextField as JTextField
import javax.swing.JButton as JButton
import javax.swing.JDialog as JDialog
import java.awt.event.ActionListener as ActionListener
import javax.swing.event.DocumentListener as DocumentListener
import javax.swing.JLabel as JLabel
import javax.swing.JPanel as JPanel
import javax.swing.SwingUtilities as SwingUtilities
import javax.swing.BorderFactory as BorderFactory
import java.awt.BorderLayout as BorderLayout
import java.awt.GridLayout as GridLayout
import java.awt.Dimension as Dimension
import javax.swing.JScrollPane as JScrollPane
import javax.swing.JSplitPane as JSplitPane
import javax.swing.JTabbedPane as JTabbedPane
import javax.swing.JTextArea as JTextArea
import javax.swing.event.ChangeListener as ChangeListener
import java.lang.Float as Float
import java.util.concurrent.ExecutionException as ExecutionException
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.VortexApplication as VortexApplication
from com.dotmatics.chemlink import ChemLink


sys.path.append(os.path.join(vortex.getVortexFolder(),"openeye"))

from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import openeye.oechem as oechem
import uk.ac.ebi.chemblws as chemblws
import uk.ac.ebi.chemblws.restclient.ChemblRestClient as ChemblRestClient
import uk.ac.ebi.chemblws.domain.Target
import org.springframework.context.support.ClassPathXmlApplicationContext as ClassPathXmlApplicationContext

CHEMBL_CMPD_LINK  = "https://www.ebi.ac.uk/chembl/compound/inspect/"
CHEMBL_ASSAY_LINK = "https://www.ebi.ac.uk/chembl/assay/inspect/"

def setSortBySelected(vortexJTable2):
    sorter = vortexJTable2.getVortexRowSorter()
    columnList = sorter.getSortedColumnList()
    orderList = sorter.getSortedOrderList()
    columnList.clear()
    columnList.add(-1) #-1 is Selection Column
    orderList.clear()
    orderList.add(0) #0 means sort by descend, 1 means sort by ascend
    sorter.fireSortingModified()
    vortexJTable2.repaint()

class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, molString, format):
        ifs = oechem.oemolistream()
        ifs.SetFormat(format)
        ifs.openstring(molString)
        mol = oechem.OEGraphMol()
        oechem.oechem.OEReadMolecule(ifs,mol)
        return mol

    def convertToInchiKey(self,mol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.INCHIKEY)
        ofs.openstring()
        oechem.oechem.OEWriteMolecule(ofs,mol)
        return ofs.GetString()

    def molToSmiles(self, mol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SMI)
        ofs.openstring()
        oechem.oechem.OEWriteMolecule(ofs,mol)
        return ofs.GetString()

class JobRunner():
    def __init__(self):
        self.progressMonitor = ProgressMonitor(vortex,"Progress","Progress", 0, 100)
        self.progressMonitor.setMillisToPopup(0)
        self.progressMonitor.setMillisToDecideToPopup(0)
        self.isCancelled = True
        self.selectedColumn = None
        self.numHits = 100
        self.nameCB = None
        self.simCutoff = 0.8
        self.useText = False
        ##################
        # from java.lang import ClassLoader
        # cl = ClassLoader.getSystemClassLoader()
        # paths = map(lambda url: url.getFile(), cl.getURLs())
        # file_io = open("/Users/jfeng1/junk.log", "w")
        # for p in paths:
        #     file_io.write(p+"\n")
        # file_io.close()
        ##################
        applicationContext = ClassPathXmlApplicationContext("applicationContext.xml")
        self.chemblClient = applicationContext.getBean("chemblRestClient")

    def show_interface(self):
        try:
            self.chemblClient.checkStatus()
        except:
            vortex.alert("Chembl webservice client is not working.")
            traceback.print_exc()
            return False
        this = self
        chemLink = ChemLink()

        if vtable is not None:
            selectedRows = vtable.getSelectedRows()

            if selectedRows.size() >= 1:
                selectedRow = selectedRows.get(0)
                mfm = vtable.getMolFileManager()
                self.queryMol = mfm.getMolFileAtRow(selectedRow)
                queryMol1 = self.queryMol
                chemLink.setMolFile(queryMol1)
                class MyRunnable(java.lang.Runnable):
                    def run(self):
                        chemLink.setMolFile(queryMol1)
                        chemLink.repaint()

                SwingUtilities.invokeLater(MyRunnable())
            chemLink.setPreferredSize(Dimension(600, 400))

        class TabChangeListener(ChangeListener):
            def stateChanged(self,e):
                tabPane = e.getSource()
                if tabPane.getSelectedIndex()==0:
                    this.useText = False
                else:
                    this.useText = True

        p = javax.swing.JPanel(BorderLayout())
        p.setBorder(BorderFactory.createTitledBorder("Click to edit structure:"))
        p.add(chemLink,BorderLayout.CENTER)

        btnPanel = JPanel()

        dbCB = JComboBox(["Chembl"])
        self.typeCB = JComboBox(["Exact","Similarity","Substructure"])
        simField = JTextField("0.8")
        simField.setEnabled(False)

        self.dialog = JDialog(VortexApplication.vortexApplication,"Database Search",True)
        self.dialog.setSize(Dimension(800,600))
        self.dialog.setLocationRelativeTo(None)
        self.dialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE)

        dialog = self.dialog
        okButton = JButton("Ok")
        class OkListener(ActionListener):
            def actionPerformed(self,e):
                this.isCancelled = False
                dialog.setVisible(False)
        okButton.addActionListener(OkListener())

        cancelButton = JButton("Cancel")
        class CancelListener(ActionListener):
            def actionPerformed(self,e):
                this.isCancelled = True
                dialog.setVisible(False)
        cancelButton.addActionListener(CancelListener())

        class SimFieldDocumentListener(DocumentListener):

            def insertUpdate(self,event):
                try:
                    this.simCutoff = float(simField.getText())
                    if this.simCutoff <=0 or this.simCutoff > 1.0:
                        okButton.setEnabled(False)
                    else:
                        okButton.setEnabled(True)
                except:
                    okButton.setEnabled(False)
                pass

            def removeUpdate(self,event):
                try:
                    this.simCutoff = float(simField.getText())
                    if this.simCutoff <=0 or this.simCutoff > 1.0:
                        okButton.setEnabled(False)
                    else:
                        okButton.setEnabled(True)
                except:
                    okButton.setEnabled(False)
                pass

            def changeUpdate(self,event):
                try:
                    this.simCutoff = float(simField.getText())
                    if this.simCutoff <=0 or this.simCutoff > 1.0:
                        okButton.setEnabled(False)
                    else:
                        okButton.setEnabled(True)
                except:
                    okButton.setEnabled(False)
                pass

        simField.getDocument().addDocumentListener(SimFieldDocumentListener())

        class TypeListener(ActionListener):
                def actionPerformed(self,e):
                    typeCB = e.getSource()
                    if typeCB.getSelectedItem()=="Similarity":
                        simField.setEnabled(True)
                    else:
                        simField.setEnabled(False)
        self.typeCB.addActionListener(TypeListener())

        tabbedPane = JTabbedPane()
        tabbedPane.addChangeListener(TabChangeListener())

        p1 = JPanel()
        p1.add(JLabel("Database:"))
        p1.add(dbCB)
        p1.add(JLabel("Type:"))
        p1.add(self.typeCB)
        p1.add(simField)
        p.add(p1,BorderLayout.SOUTH)

        tabbedPane.addTab("Structure Search",p)

        p2 = JPanel(BorderLayout())
        p2.setBorder(BorderFactory.createTitledBorder("Target/Assay Name"))
        self.area = JTextArea()
        self.area.setLineWrap(True)
        p2.add(self.area,BorderLayout.CENTER)
        tabbedPane.addTab("Assay/Target Search",p2)

        btnPanel.add(okButton)
        btnPanel.add(cancelButton)

        mainPanel = JPanel(BorderLayout())
        mainPanel.add(tabbedPane,BorderLayout.CENTER)
        mainPanel.add(btnPanel,BorderLayout.SOUTH)
        dialog.getContentPane().add(mainPanel)
        dialog.setVisible(True)
        if not this.isCancelled:
            if not self.useText:
                molfile = str(chemLink.getMolFile())
                if len(molfile.strip()) == 0:
                    return False
                else:
                    molUtilities = MolUtilities()
                    mol = molUtilities.convertToMol(molfile,oechem.OEFormat.MDL)
                    if mol.NumAtoms()>0:
                        self.inchKey = molUtilities.convertToInchiKey(mol)
                        self.smiles = molUtilities.molToSmiles(mol)
                        return True
                    else:
                        return False
            else:
                text = str(self.area.getText())
                if len(text.strip()) == 0:
                    return False
                else:
                    self.keywords = text.strip().split()
                    return True
        else:
            return False


    def run(self):
        this = self

        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self, progressMonitor):
                self.progressMonitor = progressMonitor
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)
                self.progressMonitor.setNote("Search database ...")
                self.progressMonitor.setProgress(1)
                SwingWorker.__init__(self)

            def doInBackground(self):
                compoundList = []
                searchType = "Target"
                if not this.useText:
                    searchType = this.typeCB.getSelectedItem()
                if searchType =="Exact":
                    compound = this.chemblClient.getCompoundByStdInChiKey(this.inchKey)
                    compoundList.append(compound)
                elif searchType == "Substructure":
                    cl = this.chemblClient.getCompoundBySubstructureSmiles(this.smiles.split(" ")[0])
                    print cl
                    for c in cl:
                        compoundList.append(c)
                elif searchType == "Similarity":
                    print >> sys.stderr, "Starting search"
                    cl = this.chemblClient.getSimilarCompoundBySmiles(this.smiles.split(" ")[0],int(this.simCutoff*100))
                    # print >> sys.stderr, cl
                    for c in cl:
                        compoundList.append(c)
                    print >> sys.stderr, "Done"
                else:
                    print >> sys.stderr, "getting targets ..."
                    targets = this.chemblClient.getTargets()
                    selectedTargets = []
                    for target in targets:
                        synonym = target.getSynonyms()
                        desc = target.getDescription()
                        prefName = target.getPreferredName()

                        if synonym is None:
                            synonym = ""
                        if desc is None:
                            desc = ""
                        if prefName is None:
                            prefName = ""

                        fulldesc = "%s %s %s"%(synonym,desc,prefName)
                        found = False
                        for keyword in this.keywords:
                            if fulldesc.find(keyword)>=0:
                                found = True
                                break
                        if found:
                            selectedTargets.append(target)

                if not this.useText:
                    vortexTableModel = VortexTableModel()
                    assayTableModel = VortexTableModel()
                    smiCol = vortexTableModel.findColumnWithName("Smiles",1)
                    queryCol = vortexTableModel.findColumnWithName("Query Smiles",1)
                    chemblIdCol = vortexTableModel.findColumnWithName("ChemblID",1)
                    simCol = None
                    if searchType == "Similarity":
                        simCol = vortexTableModel.findColumnWithName("Similarity",1,Column.COLUMN_DOUBLE)

                    cmpdSmilesCol_2 = assayTableModel.findColumnWithName("CmpdSmiles",1)
                    chemblIdCol_2 = assayTableModel.findColumnWithName("CmpdChemblID",1)
                    assayChemblId_2 = assayTableModel.findColumnWithName("AssayChemblID",1)
                    descriptionCol_2 = assayTableModel.findColumnWithName("Assay Description",1)
                    operatorCol_2 = assayTableModel.findColumnWithName("Operator",1)
                    valueCol_2 = assayTableModel.findColumnWithName("Value",1, Column.COLUMN_DOUBLE)
                    univCol_2 = assayTableModel.findColumnWithName("Unit",1)

                    vortexTableModel.setJoinColumn(chemblIdCol)
                    assayTableModel.setJoinColumn(chemblIdCol_2)

                    assayId = 0
                    for id,compound in enumerate(compoundList):
                        message = "Processing mol %d..."%id
                        progress = 100*id/len(compoundList)
                        self.super__publish([(message, progress)])
                        smiCol.setValueFromString(id,compound.getSmiles())
                        queryCol.setValueFromString(id,this.smiles)
                        chemblIdValue = "<html><a href=\"%s%s\">%s</a></html>" % (CHEMBL_CMPD_LINK, compound.getChemblId(), compound.getChemblId())
                        chemblIdCol.setValueFromString(id, chemblIdValue)
                        activities = this.chemblClient.getCompoundBioactivities(compound.getChemblId())
                        for act in activities:
                            cmpdSmilesCol_2.setValueFromString(assayId,compound.getSmiles())
                            chemblIdCol_2.setValueFromString(assayId,chemblIdValue)
                            assayChemblId_2.setValueFromString(assayId,"<html><a href=\"%s%s\">%s</a></html>" % (CHEMBL_ASSAY_LINK, act.getAssay_chemblid(), act.getAssay_chemblid()))
                            descriptionCol_2.setValueFromString(assayId,"<html><p>%s</p></html>"%act.getAssay_description())
                            operatorCol_2.setValueFromString(assayId,act.getOperator())
                            valueCol_2.setValueFromString(assayId,act.getValue())
                            univCol_2.setValueFromString(assayId,act.getUnits())
                            assayId += 1
                        print >> sys.stderr, activities
                        if searchType == "Similarity":
                            simCol.setValueFromString(id,"%f"%round(compound.similarity,2))
                    return vortexTableModel,assayTableModel

                else:
                    assayTableModel = VortexTableModel()
                    cmpdSmilesCol_2 = assayTableModel.findColumnWithName("CmpdSmiles",1)
                    chemblIdCol_2 = assayTableModel.findColumnWithName("CmpdChemblID",1)
                    assayChemblId_2 = assayTableModel.findColumnWithName("AssayChemblID",1)
                    descriptionCol_2 = assayTableModel.findColumnWithName("Assay Description",1)
                    operatorCol_2 = assayTableModel.findColumnWithName("Operator",1)
                    valueCol_2 = assayTableModel.findColumnWithName("Value",1, Column.COLUMN_DOUBLE)
                    univCol_2 = assayTableModel.findColumnWithName("Unit",1)
                    assayId = 0
                    for target in selectedTargets:
                        message = "Processing target %d %s..."%(assayId,target.getSynonyms())
                        progress = 100*assayId/len(selectedTargets)
                        self.super__publish([(message, progress)])
                        activities = this.chemblClient.getTargetBioactivities(target.getChemblId())
                        print >>sys.stderr, "%s %d"%(target.getSynonyms(),len(activities))
                        actId = 1
                        for act in activities:
                            message = "Processing activity %d out of %d"%(actId,len(activities))
                            self.super__publish([(message, progress)])
                            cmpd_chemblid = act.getParent_cmpd_chemblid()
                            compound = this.chemblClient.getCompound(cmpd_chemblid)
                            chemblIdValue = "<html><a href=\"%s%s\">%s</a></html>" % (CHEMBL_CMPD_LINK, compound.getChemblId(), compound.getChemblId())
                            cmpdSmilesCol_2.setValueFromString(assayId,compound.getSmiles())
                            chemblIdCol_2.setValueFromString(assayId,chemblIdValue)
                            assayChemblId_2.setValueFromString(assayId,"<html><a href=\"%s%s\">%s</a></html>" % (CHEMBL_ASSAY_LINK, act.getAssay_chemblid(), act.getAssay_chemblid()))
                            descriptionCol_2.setValueFromString(assayId,"<html><p>%s</p></html>"%act.getAssay_description())
                            operatorCol_2.setValueFromString(assayId,act.getOperator())
                            valueCol_2.setValueFromString(assayId,act.getValue())
                            univCol_2.setValueFromString(assayId,act.getUnits())
                            assayId += 1
                            actId += 1

                    return None,assayTableModel

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setNote(pair[0])
                self.progressMonitor.setProgress(pair[1])

            def done(self):
                try:
                    cmpdTblModel,assayTblModel = self.get()  #raise exception if abnormal completion
                    if cmpdTblModel is not None:
                        cmpdTbl = VortexJTable2()
                        cmpdTbl.setModel(cmpdTblModel)

                        assayTbl = VortexJTable2()
                        assayTbl.setModel(assayTblModel)
                        assayTbl.setSwitchoffColors(True)
                        setSortBySelected(assayTbl)

                        ws1 = vortex.createWorkspace(cmpdTblModel,"Chembl Compounds")
                        ws1.addComponent(JScrollPane(cmpdTbl),"Chembl Compounds")
                        ws1.updateLayout(VortexWorkspace.LAYOUT_GRID)

                        ws2 = vortex.createWorkspace(assayTblModel,"Chembl Assay Results")
                        ws2.addComponent(JScrollPane(assayTbl),"Chembl Assay Results")
                        ws2.updateLayout(VortexWorkspace.LAYOUT_GRID)

                    else:
                        ws = vortex.createWorkspace(assayTblModel,"Chembl Compounds")
                        assayTbl = VortexJTable2()
                        assayTbl.setModel(assayTblModel)
                        ws.addComponent(JScrollPane(assayTbl))
                        ws.updateLayout(VortexWorkspace.LAYOUT_GRID)

                except ExecutionException, e:
                    traceback.print_exc()
                    vortex.alert("No result found.")
                finally:
                    self.progressMonitor.close()
                return

        task = MyTask(self.progressMonitor)
        task.execute()

if vws is None:
    vortex.alert("You must have a workspace loaded...")
else:
    jobRunner = JobRunner()
    while 1:
        if jobRunner.show_interface():
            jobRunner.run()
            break
        else:
            if jobRunner.isCancelled:
                break
            else:
                if not jobRunner.useText:
                    vortex.alert("No molecule available.")
                else:
                    vortex.alert("No text available.")
