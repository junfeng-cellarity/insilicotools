package com.insilico.application.insilicotools.gui;

import chemaxon.calculations.clean.Cleaner;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import openeye.oechem.OEGraphMol;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.sql.SQLException;

/**
 * Created by jfeng1 on 1/4/16.
 */
public class CompoundInputDialog extends StandardDialog{
    MSketchPane sketcher;
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;
    Molecule templateMol;
    JTextField bioField = new JTextField(10);

    public CompoundInputDialog() {
        super();
        setModal(true);
        setSize(new Dimension(800,800));
        sketcher = MarvinFactory.createCompoundSketcher();
        okBtn = new JButton("OK");
        cancelBtn = new JButton("Cancel");
        setupTimer();
    }

    public CompoundInputDialog(JFrame parentFrame) {
        super(parentFrame);
        setModal(true);
        setSize(new Dimension(800,800));
        sketcher = MarvinFactory.createCompoundSketcher();
        okBtn = new JButton("OK");
        cancelBtn = new JButton("Cancel");
        setupTimer();
    }

    private void setupTimer(){
        Timer timer = new Timer(1000, new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(sketcher.getMol().isEmpty()){
                    okBtn.setEnabled(false);
                }else{
                    okBtn.setEnabled(true);
                }
            }
        });
        timer.start();
    }


    public boolean isCommitted() {
        return isCommitted;
    }

    public Molecule getMolecule(){
        Molecule mol = sketcher.getMol();
        if(!bioField.getText().isEmpty()){
            mol.setName(bioField.getText());
        }else{
            mol.setName(null);
        }
        return mol;
    }

    public void setChemAxonMol(Molecule mol){
        if(mol!=null){
            templateMol = mol;
            if(templateMol.getDim()!=2){
                Cleaner.clean(templateMol,2);
            }
            if(templateMol.getExplicitHcount()>0) {
                MolHandler mh = new MolHandler(templateMol);
                mh.removeHydrogens();
            }
            sketcher.setMol(templateMol);
            okBtn.setEnabled(true);
        }
    }

    public void setMolecule(PropertyMolecule mol){
        if(mol!=null){
            OEGraphMol mol1 = mol.getMol();
            templateMol = OEChemFunc.getInstance().convertOEChemMol(mol1);
            if(templateMol.getDim()!=2){
                Cleaner.clean(templateMol,2);
            }
            if(templateMol.getExplicitHcount()>0) {
                MolHandler mh = new MolHandler(templateMol);
                mh.removeHydrogens();
            }
            sketcher.setMol(templateMol);
            okBtn.setEnabled(true);
        }
    }

    @Override
    public JComponent createBannerPanel() {
        JPanel p = new JPanel();
        final JButton btn = new JButton("Load");
        btn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(bioField.getText().isEmpty()){
                    JOptionPane.showMessageDialog(CompoundInputDialog.this,"No CY Number Available.");
                    return;
                }
                String formattedCYNumber = ChemFunc.formatCYNumber(bioField.getText());
                if(formattedCYNumber==null){
                    return;
                }
                try {
                    OEGraphMol mol = InHouseCollectionDAO.getInstance().getMolFromCYNumber(formattedCYNumber);
                    if(mol!=null) {
                        bioField.setText(formattedCYNumber);
                        mol.SetTitle(formattedCYNumber);
                        sketcher.setMol(OEChemFunc.getInstance().convertOEChemMol(mol));
                        if(sketcher.getMol().isEmpty()){
                            okBtn.setEnabled(false);
                        }else{
                            okBtn.setEnabled(true);
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
        p.add(new JLabel("CY-NUMBER/Name:"));
        p.add(bioField);
        p.add(btn);
        return p;
    }

    @Override
    public JComponent createContentPanel() {
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBorder(BorderFactory.createTitledBorder("Sketch a molecule:"));
        if(templateMol!=null){
            sketcher.setMol(templateMol);
            okBtn.setEnabled(true);
        }
        sketcher.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(evt.getPropertyName().equals("mol")){
                    if(bioField.getText().startsWith("CY")) {
                        bioField.setText("");
                    }
                    if(sketcher.getMol().isEmpty()){
                        okBtn.setEnabled(false);
                    }else{
                        okBtn.setEnabled(true);
                    }
                }
            }
        });
        mainPanel.add(sketcher,BorderLayout.CENTER);
//        mainPanel.add(bottomPanel,BorderLayout.SOUTH);
        return mainPanel;
    }

    @Override
    public ButtonPanel createButtonPanel() {
        ButtonPanel buttonPanel = new ButtonPanel(SwingConstants.CENTER);
        okBtn.setEnabled(false);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(!sketcher.getMol().isEmpty()){
                    isCommitted = true;
                }else{
                    isCommitted = false; //extra layer of protection from empty molecule
                }
                CompoundInputDialog.this.setVisible(false);
            }
        });

        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                CompoundInputDialog.this.setVisible(false);
            }
        });
        buttonPanel.addButton(okBtn,ButtonPanel.AFFIRMATIVE_BUTTON);
        buttonPanel.addButton(cancelBtn,ButtonPanel.CANCEL_BUTTON);
        return buttonPanel;
    }

    public void setVisible(boolean b) {
        super.setVisible(b);
        if (b) {
            sketcher.requestFocus();
        }
    }

    public static void main(String[] args) {
        CompoundInputDialog dialog = new CompoundInputDialog();
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setVisible(true);
    }
}
