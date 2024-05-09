package com.insilico.application.insilicotools.gui.modeling;

import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.dialog.MultiInputDialog;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.gui.widget.Mol3DTablePanel;
import com.insilico.application.insilicotools.gui.widget.MolViewer3D;
import com.insilico.application.insilicotools.inSilicoTools;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.jidesoft.range.NumericRange;
import com.jidesoft.range.Range;
import openeye.oechem.OEGraphMol;

import javax.swing.*;
import javax.swing.border.BevelBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.sql.SQLException;
import java.util.*;
import java.util.List;

/**
 * Created by jfeng1 on 5/20/16.
 */
public class ConfSearchInputPanel extends JPanel{
    MSketchPane sketcher2d;
    MolViewer3D molViewer3D;
    PropertyMolecule pmol;
    JFrame resultFrame;
    Mol3DTablePanel resultPanel;
    DesignProgressMonitor progressMonitor;
    MultiInputDialog optionsDialog;
    String confMode = "fast";
    HashMap<String,Number> defaultOptions = new HashMap<String, Number>(){{
        put("No. Conformations",100);
        put("RMSD",0.1);
        put("Energy Window (kcal)", 10);
    }};
    HashMap<String,Range> defaultRanges = new HashMap<String, Range>(){{
        put("No. Conformations",new NumericRange(100,2000));
        put("RMSD",new NumericRange(0.01,2.0));
        put("Energy Window (kcal)", new NumericRange(1,20));
    }};

    public ConfSearchInputPanel() {
        super(new BorderLayout());
        progressMonitor = new DesignProgressMonitor(ConfSearchInputPanel.this,"Progress","Progress",0,100);
        optionsDialog = new MultiInputDialog((JFrame)getTopLevelAncestor(),defaultOptions,defaultRanges);
        optionsDialog.setLocationRelativeTo(this);
        optionsDialog.pack();
        JPanel p = new JPanel(new GridLayout(1,2));
        JPanel p1 = new JPanel(new BorderLayout());
        sketcher2d = MarvinFactory.createCompoundSketcher();
        sketcher2d.setSketchMode(MSketchPane.SM_SELECT_LASSO);
        sketcher2d.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(evt.getPropertyName().equals("mol")){
                    sketcher2d.centralizeMoleculeDisplay();
                    Molecule mol = sketcher2d.getMol();
                    if(!mol.isEmpty()) {
                        OEGraphMol oemol = OEChemFunc.getInstance().convertChemAxonMol(mol);
                        OEGraphMol mol3D = OEChemFunc.getInstance().getMol3D(oemol);
                        if(mol3D!=null){
                            if(pmol==null){
                                pmol = new PropertyMolecule(mol3D);
                                molViewer3D.addLigand(pmol);
                            }else{
                                pmol.setOEMol(mol3D);
                                molViewer3D.updateLigand(pmol);
                            }
                        }
                    }
                }
            }
        });
        molViewer3D = new MolViewer3D(true,false);
        p1.add(sketcher2d);

        JPanel p2 = new JPanel();
        final JTextField bioField = new JTextField(10);
        final JButton btn = new JButton("Load");
        btn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(bioField.getText().isEmpty()){
                    JOptionPane.showMessageDialog(ConfSearchInputPanel.this,"No CY Number Available.");
                    return;
                }
                String formattedCYNumber = ChemFunc.formatCYNumber(bioField.getText());
                if(formattedCYNumber==null){
                    JOptionPane.showMessageDialog(ConfSearchInputPanel.this,"Illegal CY Number Available.");
                    return;
                }
                try {
                    OEGraphMol molFromCYNumber = InHouseCollectionDAO.getInstance().getMolFromCYNumber(formattedCYNumber);
                    if(molFromCYNumber!=null) {
                        bioField.setText(formattedCYNumber);
                        Molecule mol = OEChemFunc.getInstance().convertOEChemMol(molFromCYNumber);
                        mol.setName(formattedCYNumber);
                        sketcher2d.setMol(mol);
                        btn.setEnabled(true);
                        OEGraphMol oemol = OEChemFunc.getInstance().convertChemAxonMol(mol);
                        OEGraphMol mol3D = OEChemFunc.getInstance().getMol3D(oemol);
                        if(mol3D!=null){
                            if(pmol==null){
                                pmol = new PropertyMolecule(mol3D);
                                molViewer3D.addLigand(pmol);
                            }else{
                                pmol.setOEMol(mol3D);
                                molViewer3D.updateLigand(pmol);
                            }
                        }
                    }
                } catch (SQLException e1) {
                    e1.printStackTrace();
                }
            }
        });
        btn.setEnabled(false);
        bioField.getDocument().addDocumentListener(new DocumentListener() {
            private void updateBtn(){
                if(bioField.getText().trim().isEmpty()){
                    btn.setEnabled(false);
                }else{
                    btn.setEnabled(true);
                }
            }
            @Override
            public void insertUpdate(DocumentEvent e) {
                updateBtn();
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                updateBtn();
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                updateBtn();
            }
        });
        p2.add(new JLabel("CY-NUMBER:"));
        p2.add(bioField);
        p2.add(btn);
        p1.add(p2,BorderLayout.SOUTH);

        p.add(p1);
        p.add(molViewer3D);
        add(p,BorderLayout.CENTER);
        add(buildBtnPanel(),BorderLayout.SOUTH);
        resultPanel = new Mol3DTablePanel();
        resultFrame = new JFrame("ConfSearch Result");
        resultFrame.setJMenuBar(resultPanel.getMenuBar(Mol3DTablePanel.CONF_SEARCH_MODE));
        resultFrame.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        resultFrame.setSize(new Dimension(1024,768));
        resultFrame.setLocationRelativeTo(ConfSearchInputPanel.this);
        resultFrame.getContentPane().add(resultPanel);

    }

    JPanel buildBtnPanel(){
        JPanel p = new JPanel();
        JPanel ffp = new JPanel();
        ffp.setBorder(new BevelBorder(1));
        ffp.add(new JLabel("Mode:"));
        final JComboBox ffcb = new JComboBox(new String[]{"fast","accurate"});
        ffcb.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                confMode = (String)ffcb.getSelectedItem();
            }
        });
        ffp.add(ffcb);
        p.add(ffp);
        JButton optionBtn = new JButton("Options");
        optionBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                optionsDialog.setVisible(true);
            }
        });
        p.add(optionBtn);

        JButton confSearchBtn = new JButton("Generate Conformation");
        confSearchBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                progressMonitor.setMillisToPopup(0);
                progressMonitor.setMillisToDecideToPopup(0);
                progressMonitor.setProgress(-1);
                HashMap<String,Number> map = defaultOptions;
                if(optionsDialog.isCommitted()){
                    map = optionsDialog.getResultMap();
                }
                final int numConfs = map.get("No. Conformations").intValue();
                final double rmsd = map.get("RMSD").doubleValue();
                final double ewindow = map.get("Energy Window (kcal)").doubleValue();
                inSilicoTools.getInstance().logConformationGeneration(pmol.getSmiles());
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        return OEChemFunc.getInstance().getMultiConformers(pmol.getMol3d(), numConfs, ewindow, rmsd, confMode, new ProgressReporter() {
                            @Override
                            public void reportProgress(String note, int progress) {
                                Vector v = new Vector();
                                v.add(note);
                                v.add(progress);
                                publish(v);
                            }
                        });
                    }

                    @Override
                    protected void process(List chunks) {
                        Vector v = (Vector) chunks.get(chunks.size()-1);
                        String note = (String)v.get(0);
                        int progress = (Integer) v.get(1);
                        progressMonitor.setNote(note);
                        progressMonitor.setProgress(progress);
                        if(progress==-1){
                            progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                        }
                    }

                    @Override
                    protected void done() {
                        try {
                            Vector<PropertyMolecule> result = (Vector<PropertyMolecule>) get();
                            resultPanel.setPropertyMolecules(result);
                            resultPanel.addProperty("Energy");
                            resultPanel.addProperty("rmsd_to_reference");
                            resultFrame.setVisible(true);
                        } catch (Exception e1) {
                            JOptionPane.showMessageDialog(ConfSearchInputPanel.this,e1.getMessage());
                            e1.printStackTrace();
                        }finally {
                            progressMonitor.close();
                        }
                    }
                };
                if(pmol!=null){
                    sw.execute();
                }
            }
        });
        p.add(confSearchBtn);
        return p;
    }

    public static void main(String[] args) {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(new ConfSearchInputPanel());
        frame.setSize(new Dimension(1280,1024));
        frame.setVisible(true);
    }
}
