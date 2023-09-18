package com.insilico.application.insilicotools.gui.modeling;

import chemaxon.formats.MolExporter;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.CompoundInputDialog;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.PasteSmilesDialog;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;
import com.insilico.application.insilicotools.gui.table.MolMatrixTablePane;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.github.cjwizard.StackWizardSettings;
import com.github.cjwizard.WizardContainer;
import com.github.cjwizard.WizardPage;
import com.github.cjwizard.WizardSettings;
import com.jidesoft.dialog.StandardDialog;
import openeye.oechem.*;

import javax.swing.*;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 10/18/16.
 */
public class MolWizard extends JPanel {
    WizardContainer wizard;
    WizardSettings settings;
    File currentDirectory;
    Vector<PropertyMolecule> candidates = new Vector<PropertyMolecule>();
    MatrixMolTableModel candidateTableModel = new MatrixMolTableModel(candidates);
    DesignProgressMonitor progressMonitor = new DesignProgressMonitor(MolWizard.this,"Progress","Progress",0,100);
    MolListPage molListPage = new MolListPage();
    OEGraphMol refLigand;

    public MolWizard(LayoutManager layout) {
        super(layout);
        settings = new StackWizardSettings();
    }

    private void loadSdfToTable(final Vector<PropertyMolecule> molecules, final MatrixMolTableModel tableModel){
        JFileChooser fc = new JFileChooser(currentDirectory);
        fc.setFileFilter(new FileNameExtensionFilter("SDF file.", "sdf"));
        int option = fc.showOpenDialog(MolWizard.this);
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
                            JOptionPane.showMessageDialog(MolWizard.this, e1.getMessage());
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            JOptionPane.showMessageDialog(MolWizard.this, e1.getMessage());
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

    class MolListPage extends WizardPage {
        public MolListPage() {
            super("Target", "Target Molecules");
            setLayout(new BorderLayout());
            MolMatrixTablePane tablePanel = new MolMatrixTablePane(candidateTableModel);
            add(tablePanel, BorderLayout.CENTER);
            JPanel btnPanel = new JPanel();
            JButton loadRefBtn = new JButton("Load Reference Ligand");
            loadRefBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    if(refLigand!=null){
                        candidates.insertElementAt(new PropertyMolecule(refLigand),0);
                        candidateTableModel.fireTableDataChanged();
                    }
                }
            });
            btnPanel.add(loadRefBtn);

            JButton loadSdfBtn = new JButton("Load Sdf");
            candidateTableModel.addTableModelListener(new TableModelListener() {
                @Override
                public void tableChanged(TableModelEvent e) {
                    if(candidates.size()>0){
                        MolDockingWizard.MolListPage.this.setNextEnabled(true);
                    }else{
                        MolDockingWizard.MolListPage.this.setNextEnabled(false);
                    }
                }
            });
            loadSdfBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    loadSdfToTable(candidates,candidateTableModel);
                }
            });
            btnPanel.add(loadSdfBtn);

            JButton sketchBtn = new JButton("Sketch Molecule");
            sketchBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    PropertyMolecule m = null;
                    for(PropertyMolecule mp:candidates){
                        if(mp.isSelected()){
                            m = mp;
                            break;
                        }
                    }
                    CompoundInputDialog dialog = new CompoundInputDialog();
                    if(m!=null){
                        dialog.setMolecule(m);
                    }

                    dialog.setLocationRelativeTo(MolWizard.this);
                    dialog.setVisible(true);
                    if(dialog.isCommitted()){
                        try {
                            String qmolString = MolExporter.exportToFormat(dialog.getMolecule(), "sdf");
                            OEGraphMol qmol = ChemFunc.getMolFromMolString(qmolString, OEFormat.SDF);
                            if(qmol.GetDimension()!=3) {
                                PasteSmilesDialog.processPropertyMolecules(true, candidates, qmol);
                            }else{
                                candidates.add(new PropertyMolecule(qmol));
                            }
                            candidateTableModel.fireTableDataChanged();
                        } catch (IOException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(MolWizard.this,e1.getMessage());
                        }
                    }

                }
            });
            btnPanel.add(sketchBtn);

            JButton smilesBtn = new JButton("Paste Smiles");
            smilesBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    PasteSmilesDialog dialog = new PasteSmilesDialog();
                    dialog.setLocationRelativeTo(MolWizard.this);
                    dialog.setVisible(true);
                    if(dialog.getDialogResult()== StandardDialog.RESULT_AFFIRMED){
                        Vector<PropertyMolecule> tmpMols = dialog.getMolecules();
                        if(tmpMols.size()>0){
                            if(candidates.size()>0) {
                                int option = JOptionPane.showConfirmDialog(MolWizard.this, String.format("%d molecules loaded, overwrite existing molecules?", tmpMols.size()));
                                if (option == JOptionPane.YES_OPTION) {
                                    candidates.clear();
                                }
                            }
                            for(PropertyMolecule mol:tmpMols){
                                candidates.add(mol);
                            }
                            candidateTableModel.fireTableDataChanged();
                        }
                    }
                }
            });
            btnPanel.add(smilesBtn);
            add(btnPanel,BorderLayout.SOUTH);
        }

        @Override
        public void rendering(java.util.List<WizardPage> path, WizardSettings settings) {
            super.rendering(path, settings);
            setPrevEnabled(true);
            if(candidates.size()==0) {
                setNextEnabled(false);
            }
        }
    }
}
