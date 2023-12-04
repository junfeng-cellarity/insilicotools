import sys, os
import time
import java.net.URL
import math
import jarray
import traceback

from java.lang import Runnable, Double
from java.util import Random, BitSet
from java.util.concurrent import ExecutionException

from javax.swing.border import TitledBorder
from javax.swing import BorderFactory, ProgressMonitor, JTextArea, JPanel, JLabel, JCheckBox
from javax.swing import JButton, JFrame, JTextField,JComboBox
from java.awt import GridBagLayout as awtGridBagLayout, Dimension, BorderLayout, FlowLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets
import java.awt.event.ActionListener as ActionListener
import java.awt.event.KeyListener as KeyListener
import java.awt.event.KeyEvent as KeyEvent
import java.beans.PropertyChangeListener as PropertyChangeListener
import java.lang.String as String
import com.dotmatics.vortex.table.Column as Column
from javax.swing import SwingWorker
sys.path.append(os.path.join(vortex.getVortexFolder(),"jars"))
from init_vortex import VortexInit
vortexInit = VortexInit(vortex)
vortexInit.initJars()
# vortexInit.initRDkit()
# import org.RDKit.RDKFuncs as RDKFuncs
# import org.RDKit.RWMol as RWMol
import httplib
import org.json.simple.JSONObject as JSONObject
import org.json.simple.JSONArray as JSONArray
import org.json.simple.parser.JSONParser as JSONParser
class PermPredictor:
    def __init__(self):
        self.columns = []
        for column in vtable.getColumnsArray():
            if column.isStructureColumn():
                self.columns.append(column.getName())
        self.columnCB = JComboBox(self.columns)

        if len(self.columns) == 0:
            vortex.alert("No suitable column available to convert.")
            return False

    def show_interface(self):
        mypanel = JPanel()

        columnPanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        columnLabel = JLabel("Use Column:")
        columnPanel.add(columnLabel)
        columnPanel.add(self.columnCB)

        mypanel.add(columnPanel)

        ret = vortex.showInDialog(mypanel, "Pick Structure Columnr:")
        if ret == vortex.OK:
            return True
        else:
            return False
    def calculateAll(self):
        structureCol = self.columnCB.getSelectedItem()
        this = self
        class DescriptorTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)

            def calculate_MoKa(self, smilesDict):
                conn = httplib.HTTPConnection("perm.cellarity.int", port=80)
                headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
                parser = JSONParser()

                smilesList = JSONArray()
                for key in smilesDict.keys():
                    smilesEntry = JSONObject()
                    smilesEntry.put("id",key)
                    smilesEntry.put("smiles",smilesDict[key])
                    smilesList.add(smilesEntry)
                conn.request("POST","/predict",smilesList.toJSONString(),headers)
                response = conn.getresponse()
                if response.status == 200:
                    data = response.read()
                    result_dict = parser.parse(data)
                    return result_dict
                else:
                    print(response.msg)
                    print("Something wrong.")
                return {}

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
                # {
                #     "123": {
                #         "PGP-KO_confidence": 0.007719576358795166,
                #         "in_train_PGP-KO": 0,
                #         "similarity_to_train_PGP-KO": 0.09836065769195557,
                #         "closest_smi_PGP-KO": 3.1780538303479453,
                #         "Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)": 2.094935417175293
                #     }
                # }
                batch_size = 1
                structure_column = vtable.findColumnWithName(structureCol,0)
                perm_column = vtable.findColumnWithName("Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)", vortex.DOUBLE, 1)
                smilesDict = {}
                for row in rows:
                    if self.progressMonitor.isCanceled():
                        break

                    self.super__publish([(100*row/len(rows),"process compound %d"%(row))])
                    smiles = None
                    if String(structureCol).equalsIgnoreCase("Structure") or String(structureCol).equalsIgnoreCase("MOLFILE"):
                        molfile = vtable.getMolFileManager().getMolFileAtRow(row)
                        if molfile is not None and len(molfile.strip()) > 0 and molfile.strip().upper()!="N/A":
                            try:
                                smiles = genImage.CcalcProperty(molfile,"SMILES")
                            except Exception, e:
                                smiles = None
                    else:
                        smiles = structure_column.getValueAsString(row)
                        if smiles is None or len(smiles.strip()) == 0 or smiles.strip().upper()=="N/A":
                           smiles = None

                    if smiles is not None:
                        if "\\" in smiles:
                            smiles = smiles.replace("\\","")
                        smilesDict[row] = smiles

                    if len(smilesDict)%batch_size == 0 or row == nrow -1:
                        result_dict = self.calculate_MoKa(smilesDict)
                        smilesDict.clear()
                        if len(result_dict)>0:
                            for key in result_dict.keySet():
                                rowIdx = int(key)
                                result_dict = result_dict[key]
                                perm = result_dict['Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)']
                                if perm is not None:
                                    perm_column.setValueFromString(rowIdx,"%5.2f"%float(perm))

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
        perm_runner = PermPredictor()
        if perm_runner.show_interface():
            perm_runner.calculateAll()