package com.insilico.application.insilicotools.gui.modeling;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.CompoundInputDialog;
import com.insilico.application.insilicotools.gui.dialog.MultiInputDialog;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;
import com.insilico.application.insilicotools.gui.table.MolMatrixTablePane;
import com.insilico.application.insilicotools.gui.widget.Mol3DTablePanel;
import com.insilico.application.insilicotools.gui.widget.MolViewer3D;
import com.insilico.application.insilicotools.inSilicoTools;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.util.SuperpositionSolution;
import com.github.cjwizard.*;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.jidesoft.range.NumericRange;
import com.jidesoft.range.Range;
import openeye.oechem.*;
import org.jdesktop.swingx.JXButton;
import org.jdesktop.swingx.painter.RectanglePainter;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 6/23/16.
 */
public class MolOverlayWizard extends MolWizard {
    MolViewer3D molViewer3D = new MolViewer3D(true,false);
    MolViewer3D molViewer3D_2 = new MolViewer3D(false,false);

    MatrixMolTableModel previewTableModel = new MatrixMolTableModel(candidates);

    TemplatePage templatePage = new TemplatePage();
    OverlayPage overlayPage = new OverlayPage();

    OverlayPageFactory pageFactory;

    JFrame resultFrame;
    Mol3DTablePanel resultPanel;
    CompoundInputDialog structureInputDialog;
    MultiInputDialog optionsDialog;
    HashMap<String,Number> defaultOptions = new HashMap<String, Number>(){{
        put("No. Conformations",500);
        put("RMSD",0.1);
        put("Energy Window (kcal)", 10);
    }};
    HashMap<String,Range> defaultRanges = new HashMap<String, Range>(){{
        put("No. Conformations",new NumericRange(100,10000));
        put("RMSD",new NumericRange(0.01,2.0));
        put("Energy Window (kcal)", new NumericRange(1,20));
    }};
    OverlayInputDialog overlayInputDialog;

    public MolOverlayWizard() {
        super(new BorderLayout());
        optionsDialog = new MultiInputDialog((JFrame)getTopLevelAncestor(),defaultOptions,defaultRanges);
        optionsDialog.setLocationRelativeTo(this);
        optionsDialog.pack();
        structureInputDialog = new CompoundInputDialog((JFrame)this.getTopLevelAncestor());
        structureInputDialog.setLocationRelativeTo(this);
        settings = new StackWizardSettings();
        pageFactory = new OverlayPageFactory();
        wizard = new WizardContainer(pageFactory){
            @Override
            public void finish() {
                super.finish();
                MolOverlayWizard.this.firePropertyChange("Finished",true,false);
            }

            @Override
            public void cancel() {
                super.cancel();
                MolOverlayWizard.this.firePropertyChange("Cancelled",true,false);
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
        resultPanel = new Mol3DTablePanel();
        resultFrame = new JFrame("Structure Overlay Result");
        resultFrame.setJMenuBar(resultPanel.getMenuBar(Mol3DTablePanel.SUPERIMPOSE_MODE));
        resultFrame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        resultFrame.setSize(new Dimension(1024,768));
        resultFrame.setLocationRelativeTo(MolOverlayWizard.this);
        resultFrame.getContentPane().add(resultPanel);
    }

    class TemplatePage extends WizardPage {
        JPanel btnPanel;
        public TemplatePage() {
            super("Template", "Template molecule");
            setLayout(new BorderLayout());
            add(molViewer3D,BorderLayout.CENTER);
            btnPanel = new JPanel();
            JButton loadBtn = new JButton("Load Template");
            loadBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent actionEvent) {
                    try {
                        if(overlayInputDialog==null) {
                            Vector<String> templateSructures = ChemFunc.getTemplateStructures();
                            overlayInputDialog = new OverlayInputDialog((JFrame) MolOverlayWizard.this.getTopLevelAncestor(),templateSructures);
                            overlayInputDialog.setLocationRelativeTo(MolOverlayWizard.this);
                        }
                        overlayInputDialog.setVisible(true);
                        if(overlayInputDialog.isCommitted()){
                            molViewer3D.clear();
                            molViewer3D_2.clear();
                            String ligandMolStr = overlayInputDialog.getLigandString();
                            refLigand = ChemFunc.getMolFromMolString(ligandMolStr,OEFormat.MDL);
                            molViewer3D.addLigand(new PropertyMolecule(refLigand));
                            molViewer3D_2.addLigand(new PropertyMolecule(refLigand));
                            setNextEnabled(true);
                        }
                    } catch (Exception e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MolOverlayWizard.this,e1.getMessage());
                    }


                }
            });
            btnPanel.add(loadBtn);
            JButton pasteBtn = new JButton("Paste Structure");
            pasteBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                    try {
                        Object data = clipboard.getData(DataFlavor.stringFlavor);
                        if(data!=null){
                            String sdfString = (String)data;
                            OEGraphMol mol = ChemFunc.getMolFromMolString(sdfString, OEFormat.MDL);
                            if(mol!=null) {
                                if (mol.GetDimension() != 3) {
                                    refLigand = OEChemFunc.getInstance().getMol3D(mol);
                                }else{
                                    refLigand = mol;
                                }
                            }
                            if(refLigand!=null){
                                molViewer3D.clear();
                                PropertyMolecule propertyMolecule = new PropertyMolecule(refLigand);
                                molViewer3D.addLigand(propertyMolecule);
                                molViewer3D_2.clear();
                                molViewer3D_2.addLigand(propertyMolecule);
                                TemplatePage.this.setNextEnabled(true);
                            }else{
                                TemplatePage.this.setNextEnabled(false);
                            }

                        }
                    } catch (UnsupportedFlavorException e1) {
                        e1.printStackTrace();
                    } catch (IOException e1) {
                        e1.printStackTrace();
                    }
                }
            });
            btnPanel.add(pasteBtn);

            JButton sketchBtn = new JButton("Sketch");
            sketchBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    structureInputDialog.setVisible(true);
                    if(structureInputDialog.isCommitted()&&structureInputDialog.getMolecule()!=null){
                        OEGraphMol mol = OEChemFunc.getInstance().convertChemAxonMol(structureInputDialog.getMolecule());
                        if(mol!=null) {
                            if (mol.GetDimension() != 3) {
                                refLigand = OEChemFunc.getInstance().getMol3D(mol);
                            }else{
                                refLigand = mol;
                            }
                        }
                        if(refLigand!=null){
                            molViewer3D.clear();
                            PropertyMolecule propertyMolecule = new PropertyMolecule(refLigand);
                            molViewer3D.addLigand(propertyMolecule);
                            molViewer3D_2.clear();
                            molViewer3D_2.addLigand(propertyMolecule);
                            TemplatePage.this.setNextEnabled(true);
                        }else{
                            TemplatePage.this.setNextEnabled(false);
                        }

                    }
                }
            });
            btnPanel.add(sketchBtn);
            add(btnPanel,BorderLayout.SOUTH);
        }

        @Override
        public void rendering(List<WizardPage> path, WizardSettings settings) {
            super.rendering(path, settings);
            setPrevEnabled(false);
            if(refLigand==null) {
                setNextEnabled(false);
            }
        }
    }

    class OverlayPage extends WizardPage{
        public OverlayPage() {
            super("Overlap","Overlay Molecules on Template Molecule");
            setLayout(new BorderLayout());
            JPanel p = new JPanel(new GridLayout(1,2));
            MolMatrixTablePane tablePanel = new MolMatrixTablePane(previewTableModel,true);
            p.add(molViewer3D_2);
            p.add(tablePanel);
            JPanel btnPanel = new JPanel();
            //final JCheckBox batchCb = new JCheckBox("Batch",false);
            //final JComboBox methodCb = new JComboBox(new String[]{"ROCS","EON"});
            final JButton optionBtn = new JButton("ConfGen Options");
            final JComboBox nSolutionCB = new JComboBox(new Integer[]{1,2,3,4,5,6,7,8,9,10});
            //HashMap<String,String> methodHash = new HashMap<String,String>(){{put("ROCS","oe_rocs");put("EON","oe_eon");}};
            //methodCb.setEnabled(false);
            //btnPanel.add(batchCb);
            //btnPanel.add(methodCb);
            final JCheckBox rigidCb = new JCheckBox("Rigid Overlay?",false);
            rigidCb.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent actionEvent) {
                    if(rigidCb.isSelected()) {
//                        batchCb.setEnabled(false);
//                        methodCb.setEnabled(false);
                        optionBtn.setEnabled(false);
                        nSolutionCB.setEnabled(false);
                    }else{
//                        batchCb.setEnabled(true);
//                        methodCb.setEnabled(true);
                        optionBtn.setEnabled(true);
                        nSolutionCB.setEnabled(true);
                    }
                }
            });
            btnPanel.add(rigidCb);
            //final JCheckBox shapeOnlyCb = new JCheckBox("Shape Only?",false);
            //btnPanel.add(shapeOnlyCb);
            btnPanel.add(new JLabel("No. Solution:"));
            btnPanel.add(nSolutionCB);
            optionBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    optionsDialog.setVisible(true);
                }
            });
            btnPanel.add(optionBtn);
            /*
            batchCb.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent actionEvent) {
                    if(batchCb.isSelected()){
                        nSolutionCB.setEnabled(false);
                        methodCb.setEnabled(true);
                        rigidCb.setEnabled(false);
                    }else{
                        nSolutionCB.setEnabled(true);
                        methodCb.setEnabled(false);
                        rigidCb.setEnabled(true);
                    }
                }
            });
            */
            JXButton superimposeBtn = new JXButton("Superimpose");
            superimposeBtn.setBackgroundPainter(new RectanglePainter(Color.GREEN,Color.BLACK));

            superimposeBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    progressMonitor.setMillisToDecideToPopup(1);
                    progressMonitor.setMillisToPopup(1);
                    progressMonitor.setNote("Superimposing molecules ...");
                    progressMonitor.setProgress(1);
                    inSilicoTools.getInstance().logSuperposition(oechem.OEMolToSmiles(refLigand),candidates.size());
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            NumberFormat nf = new DecimalFormat("#.##");
                            Vector<PropertyMolecule> results = new Vector<PropertyMolecule>();
                            results.add(new PropertyMolecule(refLigand));
                            if(rigidCb.isSelected()){
                                results.addAll(ChemFunc.rigid_overlay(refLigand, candidates, new ProgressReporter() {
                                    @Override
                                    public void reportProgress(String note, int progress) {
                                        Vector v = new Vector();
                                        v.add(note);
                                        v.add(progress);
                                        publish(v);
                                    }
                                }));

//                                for(PropertyMolecule candidate:candidates) {
//                                    OEGraphMol result = OEChemFunc.getInstance().colorOverlapRigid2Rigid(refLigand,candidate.getMol3d(),shapeOnlyCb.isSelected());
//                                    candidate.setConf3d(result);
//                                    if(oechem.OEHasSDData(result,"Shape/Feature Score")) {
//                                        candidate.addProperty("Overlay Score", nf.format(Double.parseDouble(oechem.OEGetSDData(result, "Shape/Feature Score"))));
//                                    }
//                                    if(oechem.OEHasSDData(result,"Energy")) {
//                                        candidate.addProperty("Strain Energy", nf.format(Double.parseDouble(oechem.OEGetSDData(result, "Energy"))));
//                                    }
//                                    results.add(candidate);
//                                }
                            }else /*if(batchCb.isSelected())*/{
                                //final String method = methodHash.get((String)methodCb.getSelectedItem());
                                HashMap<String, Number> optionMap = defaultOptions;
                                if(optionsDialog.isCommitted()){
                                    optionMap = optionsDialog.getResultMap();
                                }
                                final double rmsd = optionMap.get("RMSD").doubleValue();
                                final double ewindow = optionMap.get("Energy Window (kcal)").doubleValue();
                                final int nSolutions = (int) nSolutionCB.getSelectedItem();

                                results.addAll(ChemFunc.batch_overlay(refLigand, candidates, ewindow, rmsd, nSolutions, new ProgressReporter() {
                                    @Override
                                    public void reportProgress(String note, int progress) {
                                        Vector v = new Vector();
                                        v.add(note);
                                        v.add(progress);
                                        publish(v);
                                    }
                                }));
                            }
                            /*
                            else {
                                int nSolutions = (Integer) nSolutionCB.getSelectedItem();
                                HashMap<String, Number> optionMap = defaultOptions;
                                if(optionsDialog.isCommitted()){
                                    optionMap = optionsDialog.getResultMap();
                                }
                                final int numConfs = optionMap.get("No. Conformations").intValue();
                                final double rmsd = optionMap.get("RMSD").doubleValue();
                                final double ewindow = optionMap.get("Energy Window (kcal)").doubleValue();

                                if (candidates.size() == 1) {
                                    Vector<SuperpositionSolution> solutions = OEChemFunc.getInstance().colorOverlapRigid2Flex(refLigand, candidates.get(0).getMol3d(), numConfs, rmsd, ewindow, new ProgressReporter() {
                                        @Override
                                        public void reportProgress(String note, int progress) {
                                            Vector v = new Vector();
                                            v.add(note);
                                            v.add(progress);
                                            publish(v);
                                        }
                                    },nSolutions,shapeOnlyCb.isSelected());

                                    for (SuperpositionSolution solution : solutions) {
                                        PropertyMolecule m1 = new PropertyMolecule(solution.getTarget());
                                        m1.addProperty("Strain Energy", nf.format(solution.getEnergy()));
                                        m1.addProperty("Overlay Score", nf.format(solution.getScore()));
                                        results.add(m1);
                                    }
                                } else {
                                    int n = 0;
                                    for (PropertyMolecule m : candidates) {
                                        Vector<SuperpositionSolution> solutions = OEChemFunc.getInstance().colorOverlapRigid2Flex(
                                                refLigand, m.getMol3d(), numConfs, rmsd, ewindow,null, nSolutions, shapeOnlyCb.isSelected());
                                        int progress = 100 * n / candidates.size();
                                        publish(progress);
                                        for (SuperpositionSolution solution : solutions) {
                                            PropertyMolecule m1 = new PropertyMolecule(solution.getTarget());
                                            m1.addProperty("Strain Energy", nf.format(solution.getEnergy()));
                                            m1.addProperty("Overlay Score", nf.format(solution.getScore()));
                                            results.add(m1);
                                        }
                                        n += 1;
                                    }
                                }
                            }
                             */
                            return results;
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
                                progressMonitor.setNote("Superimposing ...");
                                progressMonitor.setProgress(progress);
                            }
                        }

                        @Override
                        protected void done() {
                            try {
                                Object result = get();
                                if(result!=null){
                                    Vector<PropertyMolecule> resultMols = (Vector<PropertyMolecule>)result;
                                    resultPanel.setPropertyMolecules(resultMols);
                                    resultPanel.addProperty("Strain Energy");
                                    resultPanel.addProperty("Overlay Score");
                                    resultFrame.setVisible(true);
                                }
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                            } finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    sw.execute();
                    //OEChemFunc.getInstance().colorOverlapRigid2Flex()
                }
            });
            btnPanel.add(superimposeBtn);
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

    private class OverlayPageFactory implements PageFactory {
        public OverlayPageFactory() {

        }

        @Override
        public WizardPage createPage(List<WizardPage> list, WizardSettings wizardSettings) {
            switch (list.size()){
                case 0:
                    return templatePage;
                case 1:
                    return molListPage;
                case 2:
                    return overlayPage;
                default:
                    return null;
            }
        }

        @Override
        public boolean isTransient(List<WizardPage> list, WizardSettings wizardSettings) {
            return false;
        }

    }

    private void loadSdfToTable(final Vector<PropertyMolecule> molecules, final MatrixMolTableModel tableModel){
        JFileChooser fc = new JFileChooser(currentDirectory);
        fc.setFileFilter(new FileNameExtensionFilter("SDF file.", "sdf"));
        int option = fc.showOpenDialog(MolOverlayWizard.this);
        if (option == JFileChooser.APPROVE_OPTION) {
            currentDirectory = fc.getCurrentDirectory();
            File f = fc.getSelectedFile();
            if (f.canRead()) {
                molecules.clear();
                final String fname = f.getAbsolutePath();
//                    OEMolDatabase db = new OEMolDatabase(fname);
//                    final int numMols = db.GetMaxMolIdx();
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        oemolistream ifs = new oemolistream();
                        ifs.SetFormat(OEFormat.SDF);
                        ifs.open(fname);
                        Vector<String> properties = new Vector<String>();
                        OEGraphMol mol = new OEGraphMol();
                        while (oechem.OEReadMolecule(ifs, mol)) {
                            OESDDataIter iter = oechem.OEGetSDDataPairs(mol);
                            while(iter.hasNext()){
                                OESDDataPair pair = iter.next();
                                if(!properties.contains(pair.GetTag())) {
                                    properties.add(pair.GetTag());
                                }
                            }
                            molecules.add(new PropertyMolecule(mol));
                        }
//                        ChemFunc.calculateOEProperty(molecules);
//                        ChemFunc.generateCLogP(molecules);
//
                        return properties;
                    }

                    @Override
                    protected void done() {
                        try {
                            Vector<String> properties = (Vector<String>)get();
//                            String[] tags = properties.toArray(new String[properties.size()]);
//                            if(tags.length>0) {
//                                String selectedTag = (String) JOptionPane.showInputDialog(MolOverlayWizard.this, "Input", "Mol Name", JOptionPane.QUESTION_MESSAGE, null, tags, tags[0]);
//                                for (PropertyMolecule m : molecules) {
//                                    String title = oechem.OEGetSDData(m.getMol(), selectedTag);
//                                    m.setName(title);
//                                }
//                            }
                            tableModel.fireTableStructureChanged();
                        } catch (InterruptedException e1) {
                            JOptionPane.showMessageDialog(MolOverlayWizard.this, e1.getMessage());
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            JOptionPane.showMessageDialog(MolOverlayWizard.this, e1.getMessage());
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();

            }
        }

    }



    public static void main(String[] args) {
        MolOverlayWizard wizard = new MolOverlayWizard();
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(wizard);
        f.setSize(new Dimension(1024,768));
        f.setVisible(true);
    }

}
