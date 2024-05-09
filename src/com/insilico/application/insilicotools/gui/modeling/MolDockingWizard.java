package com.insilico.application.insilicotools.gui.modeling;

import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;
import com.insilico.application.insilicotools.gui.table.MolMatrixTablePane;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.gui.widget.Mol3DTablePanel;
import com.insilico.application.insilicotools.gui.widget.MolViewer3D;
import com.insilico.application.insilicotools.inSilicoTools;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.github.cjwizard.*;
import com.insilico.application.insilicotools.gui.CompoundInputDialog;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.SubstructureSelectionDialog;
import openeye.oechem.*;
import org.jdesktop.swingx.JXButton;
import org.jdesktop.swingx.painter.RectanglePainter;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import openeye.oedocking.*;

/**
 * Created by jfeng1 on 6/23/16.
 */
public class MolDockingWizard extends MolWizard {
    MolViewer3D molViewer3D = new MolViewer3D(true,true);
    MolViewer3D molViewer3D_2 = new MolViewer3D(false,false);

    MatrixMolTableModel previewTableModel = new MatrixMolTableModel(candidates);

    ReceptorPage receptorPage = new ReceptorPage();
    DockingPage dockingPage = new DockingPage();

    String receptorPdbStr;
    String ligandMolStr;
    String receptorName;
    DockingPageFactory pageFactory;

    JFrame resultFrame;
    Mol3DTablePanel resultPanel;
    CompoundInputDialog structureInputDialog;
    DockingInputDialog dockingInputDialog;
    String atomIdStr;
    String smarts;
    Molecule refMol;
    String refLigandStr;
    String core_smarts;
    String core_atoms;
    String allowed_clashes;
//          -allowed_clashes : Allowed severity of clashes, values may be one of
//                         "noclashes" ( less than 0.2A penetration), "mildclashes" ( less
//    than 0.65A penetration), and "allclashes"


    public MolDockingWizard() {
        super(new BorderLayout());
        structureInputDialog = new CompoundInputDialog((JFrame)this.getTopLevelAncestor());
        structureInputDialog.setLocationRelativeTo(this);
        pageFactory = new DockingPageFactory();
        wizard = new WizardContainer(pageFactory){
            @Override
            public void finish() {
                super.finish();
                MolDockingWizard.this.firePropertyChange("Finished",true,false);
            }

            @Override
            public void cancel() {
                super.cancel();
                MolDockingWizard.this.firePropertyChange("Cancelled",true,false);
            }
        };
        wizard.addWizardListener(new WizardListener() {
            @Override
            public void onPageChanged(WizardPage wizardPage, List<WizardPage> list) {

            }

            @Override
            public void onFinished(List<WizardPage> list, WizardSettings wizardSettings) {

            }

            @Override
            public void onCanceled(List<WizardPage> list, WizardSettings wizardSettings) {

            }
        });
        add(wizard,BorderLayout.CENTER);
        resultPanel = new Mol3DTablePanel(true);
        resultFrame = new JFrame("Docking Result");
        resultFrame.setJMenuBar(resultPanel.getMenuBar(Mol3DTablePanel.SUPERIMPOSE_MODE));
        resultFrame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        resultFrame.setSize(new Dimension(1024,768));
        resultFrame.setLocationRelativeTo(MolDockingWizard.this);
        resultFrame.getContentPane().add(resultPanel);
    }

    class ReceptorPage extends WizardPage {
        JPanel btnPanel;
        public ReceptorPage() {
            super("Docking Structure", "Docking Structure");
            setLayout(new BorderLayout());
            add(molViewer3D,BorderLayout.CENTER);
            btnPanel = new JPanel();
            JButton loadBtn = new JButton("Load Structure");
            loadBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    try {
                        if(dockingInputDialog==null) {
                            Vector<String> dockingStructures = ChemFunc.getDockingStructures();
                            dockingInputDialog = new DockingInputDialog((JFrame) MolDockingWizard.this.getTopLevelAncestor(), dockingStructures);
                            dockingInputDialog.setLocationRelativeTo(MolDockingWizard.this);
                        }
                        dockingInputDialog.setVisible(true);
                        if(dockingInputDialog.isCommitted()){
                            molViewer3D.clear();
                            receptorName = dockingInputDialog.getReceptorName();
                            receptorPdbStr = dockingInputDialog.getPdbString();
                            ligandMolStr = dockingInputDialog.getLigandString();
                            refLigand = ChemFunc.getMolFromMolString(ligandMolStr,OEFormat.MDL);
                            refMol = null;
                            refLigandStr = null;
                            core_atoms = null;
                            core_smarts = null;
                            molViewer3D.setReceptor(receptorName,receptorPdbStr);
                            molViewer3D.addLigand(new PropertyMolecule(refLigand));
                            molViewer3D_2.setReceptor(receptorName,receptorPdbStr);
                            molViewer3D_2.addLigand(new PropertyMolecule(refLigand));
                            setNextEnabled(true);
                        }
                    } catch (Exception e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MolDockingWizard.this,e1.getMessage());
                    }
                }
            });
            btnPanel.add(loadBtn);

            add(btnPanel,BorderLayout.SOUTH);
        }

        @Override
        public void rendering(List<WizardPage> path, WizardSettings settings) {
            super.rendering(path, settings);
            setPrevEnabled(false);
            if(ligandMolStr==null||receptorPdbStr==null) {
                setNextEnabled(false);
            }
        }
    }

    class DockingPage extends WizardPage{
        public DockingPage() {
            super("Docking","Docking Molecules on selected receptor");
            setLayout(new BorderLayout());
            JPanel p = new JPanel(new GridLayout(1,2));
            final MolMatrixTablePane tablePanel = new MolMatrixTablePane(previewTableModel,true);
            p.add(molViewer3D_2);
            p.add(tablePanel);
            JPanel btnPanel = new JPanel();
//            final JComboBox dockingModeCB = new JComboBox(new String[]{"Glide:SP:Faster","Glide:XP:Extensive","Glide:Rigid","OpenEye:hybrid","Openeye:posit"});
            final JComboBox dockingModeCB = new JComboBox(new String[]{"Glide:SP:Faster","Glide:XP:Extensive","Glide:Rigid"});
//            final JComboBox clashesCB = new JComboBox(new String[]{"allclashes","mildclashes","noclashes"});
//            clashesCB.setEnabled(false);


            btnPanel.add(new JLabel("Mode"));
            btnPanel.add(dockingModeCB);
            final JComboBox numPosesCB = new JComboBox(new Integer[]{1,2,3,4,5,6,7,8,9,10});
            btnPanel.add(new JLabel("No. Poses:"));
            btnPanel.add(numPosesCB);
            //btnPanel.add(clashesCB);

            final JCheckBox consCB = new JCheckBox("Use Constraint",false);
            final JCheckBox macrocycleCB = new JCheckBox("All Macrocycles?", false);
            final JButton constraintBtn = new JButton("Define Constraint");
            constraintBtn.setEnabled(false);

            dockingModeCB.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    String selected = (String)dockingModeCB.getSelectedItem();
                    String vendor = selected.split(":")[0];
                    String docking_method = selected.split(":")[1];
                    if(vendor.equals("Glide")){
                        consCB.setEnabled(true);
                        if(consCB.isSelected()){
                            constraintBtn.setEnabled(true);
                        }
                        //clashesCB.setEnabled(false);
                    }else{
                        consCB.setEnabled(false);
                        constraintBtn.setEnabled(false);
//                        if(docking_method.equalsIgnoreCase("posit")){
//                            clashesCB.setEnabled(true);
//                        }else{
//                            clashesCB.setEnabled(false);
//                        }
                    }
                }
            });

            final JCheckBox protonationCB = new JCheckBox("Protonation by pKa?", true);

            JXButton dockingBtn = new JXButton("Dock");
            dockingBtn.setBackgroundPainter(new RectanglePainter(Color.GREEN,Color.RED));
            dockingBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    progressMonitor.setMillisToDecideToPopup(1);
                    progressMonitor.setMillisToPopup(1);
                    progressMonitor.setNote("Docking molecules ...");
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    final int numPoses = (Integer)numPosesCB.getSelectedItem();
                    String[] args = ((String) dockingModeCB.getSelectedItem()).split(":");
                    final String dockingVendor = args[0];
                    final String dockingMode;
                    final String dockingMethod;
                    if(dockingVendor.equals("Glide")) {
                        if (args[1].equals("Rigid")) {
                            dockingMode = "SP";
                            dockingMethod = "rigid";
                        } else {
                            dockingMode = args[1];
                            dockingMethod = "confgen";
                        }
                    }else{
                        dockingMode = "";
                        dockingMethod = args[1];
                    }
                    inSilicoTools.getInstance().logDocking(receptorName,candidates.size());
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            if(dockingVendor.equals("Glide")) {
                                return ChemFunc.dock(receptorName, candidates, numPoses, dockingMethod, dockingMode, macrocycleCB.isSelected(),
                                        consCB.isSelected(), refLigandStr == null ? ligandMolStr : refLigandStr, core_atoms, core_smarts, protonationCB.isSelected(),
                                        new ProgressReporter() {
                                            @Override
                                            public void reportProgress(String note, int progress) {
                                                Vector v = new Vector();
                                                v.add(note);
                                                v.add(progress);
                                                publish(v);
                                            }
                                        });
                            }
//                            else {
//                                return ChemFunc.oedock(receptorName, candidates, dockingMethod.equals("posit")?dockingMethod+"_"+clashesCB.getSelectedItem():dockingMethod,
//                                        protonationCB.isSelected(), refLigandStr == null ? ligandMolStr : refLigandStr, numPoses,
//                                        new ProgressReporter() {
//                                            @Override
//                                            public void reportProgress(String note, int progress) {
//                                                Vector v = new Vector();
//                                                v.add(note);
//                                                v.add(progress);
//                                                publish(v);
//                                            }
//                                        });
//                            }
                            return null;
                        }

                        @Override
                        protected void process(List chunks) {
                            Object obj = chunks.get(chunks.size()-1);
                            if(obj instanceof Vector){
                                Vector v = (Vector)obj;
                                String note = (String)v.get(0);
                                int progress = (Integer)v.get(1);
                                progressMonitor.setNote(note);
                                progressMonitor.setProgress(progress);
                            }else{
                                int progress = (Integer)obj;
                                progressMonitor.setNote("Docking molecules ...");
                                progressMonitor.setProgress(progress);
                            }
                        }

                        @Override
                        protected void done() {
                            try {
                                Object result = get();
                                if(result!=null){
                                    Vector<PropertyMolecule> resultMols = (Vector<PropertyMolecule>)result;
                                    PropertyMolecule ref = new PropertyMolecule(ChemFunc.getMolFromMolString(ligandMolStr,OEFormat.MDL));
                                    resultMols.insertElementAt(ref,0);
                                    resultPanel.setPropertyMolecules(receptorName,receptorPdbStr,resultMols);
                                    if(dockingVendor.equals("Glide")) {
                                        resultPanel.removeProperty("Probability");
                                        resultPanel.removeProperty("ChemGauss4 Score");
                                        resultPanel.addProperty("Docking Score");
                                        resultPanel.addProperty("Overlay Score");
                                        resultPanel.addProperty("Interaction Energy");
                                    }else{
                                        resultPanel.removeProperty("Docking Score");
                                        resultPanel.removeProperty("Interaction Energy");
                                        if(dockingMethod.equalsIgnoreCase("posit")){
                                            resultPanel.removeProperty("ChemGauss4 Score");
                                            resultPanel.addProperty("Overlay Score");
                                            resultPanel.addProperty("Probability");
                                        }else {
                                            resultPanel.removeProperty("Probability");
                                            resultPanel.addProperty("ChemGauss4 Score");
                                            resultPanel.addProperty("Overlay Score");
                                        }
                                    }
//                                    resultPanel.addProperty("E_Internal");
                                    progressMonitor.close();
                                    resultFrame.setVisible(true);
                                }
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                                progressMonitor.close();
                                JOptionPane.showMessageDialog(MolDockingWizard.this,e1.getMessage());
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(MolDockingWizard.this,e1.getMessage());
                                progressMonitor.close();
                            } finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    sw.execute();
                }
            });
            btnPanel.add(dockingBtn);

            btnPanel.add(protonationCB);

            btnPanel.add(consCB);
            btnPanel.add(macrocycleCB);
            constraintBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    SubstructureSelectionDialog dialog = new SubstructureSelectionDialog();
                    dialog.setLocationRelativeTo(MolDockingWizard.this);
                    if(refMol==null) {
                        dialog.setPropertyMolecule(new PropertyMolecule(refLigand));
                    }else{
                        dialog.setMolecule(refMol);
                    }
                    dialog.setVisible(true);
                    if(!dialog.isCommitted()){
                        return;
                    }
                    refMol = dialog.getMolecule();
                    core_smarts = dialog.getSmarts();
                    core_atoms = dialog.getAtomIdsStr();
                    oechem.OEAddExplicitHydrogens(refLigand);
                    refLigandStr = ChemFunc.getMolString(refLigand);
                }
            });
            btnPanel.add(constraintBtn);
            consCB.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    constraintBtn.setEnabled(consCB.isSelected());
                }
            });
            add(p,BorderLayout.CENTER);
            add(btnPanel,BorderLayout.SOUTH);
        }

        @Override
        public void rendering(List<WizardPage> path, WizardSettings settings) {
            super.rendering(path, settings);
            setPrevEnabled(true);
            setNextEnabled(false);
            setFinishEnabled(false);
            previewTableModel.fireTableStructureChanged();
            previewTableModel.fireTableDataChanged();
        }
    }

    private class DockingPageFactory implements PageFactory {
        public DockingPageFactory() {

        }

        @Override
        public WizardPage createPage(List<WizardPage> list, WizardSettings wizardSettings) {
            atomIdStr = null;
            smarts = null;
            switch (list.size()){
                case 0:
                    return receptorPage;
                case 1:
                    return molListPage;
                case 2:
                    return dockingPage;
                default:
                    return null;
            }
        }

        @Override
        public boolean isTransient(List<WizardPage> list, WizardSettings wizardSettings) {
            return false;
        }

    }


    public static void main(String[] args) {
        JPopupMenu.setDefaultLightWeightPopupEnabled(false);
        try {
            File f = new File("/Users/jfeng1/oe_license.txt");
            oechem.OEAddLicenseData(FileUtil.readFileToString(f));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        String receptor_file = "/Users/jfeng1/Datasets/PostIt/receptor.pdb";
        String ref_file = "/Users/jfeng1/Datasets/PostIt/ligand.sdf";
        String lig_file = "/Users/jfeng1/Datasets/PostIt/unknown3.sdf";
        oemolistream ifs = new oemolistream();
        ifs.open(ref_file);
        OEGraphMol ref_mol = new OEGraphMol();
        oechem.OEReadMolecule(ifs,ref_mol);
        ifs.close();

        ifs.open(lig_file);
        OEGraphMol ligand_mol = new OEGraphMol();
        oechem.OEReadMolecule(ifs,ligand_mol);
        ifs.close();

        OEMol mcMol = new OEMol(ligand_mol);
        OEChemFunc.getInstance().generateMultipleConformers(mcMol,1000,10);

        ifs.open(receptor_file);
        OEGraphMol protein = new OEGraphMol();
        oechem.OEReadMolecule(ifs,protein);
        ifs.close();
        OEGraphMol receptor = new OEGraphMol();
        oedocking.OEMakeReceptor(receptor, protein, ref_mol);
        OEPositOptions option = new OEPositOptions();
        option.SetFullConformationSearch(true);
        option.SetPositMethods(OEPositMethod.SHAPEFIT);
        OEPosit dock = new OEPosit(option);
//        OEDock dock = new OEDock(OEDockMethod.Hybrid2, OESearchResolution.High);
        dock.Initialize(receptor);
        OEGraphMol dockedMol = new OEGraphMol();
        dock.DockMultiConformerMolecule(dockedMol, mcMol);
        oemolostream ofs = new oemolostream();
        ofs.open("/Users/jfeng1/Datasets/PostIt/docked.sdf");
        dock.DockMultiConformerMolecule(dockedMol,mcMol);
        dock.ScoreLigand(dockedMol);
        System.out.println(dockedMol.GetEnergy());
        oechem.OEWriteMolecule(ofs,dockedMol);
        ofs.close();

        /*
        MolDockingWizard wizard = new MolDockingWizard();
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(wizard);
        f.setSize(new Dimension(1024,768));
        f.setVisible(true);
        */
    }

}
