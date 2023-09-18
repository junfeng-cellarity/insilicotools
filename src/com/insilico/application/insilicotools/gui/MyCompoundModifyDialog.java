package com.insilico.application.insilicotools.gui;

import chemaxon.formats.MolExporter;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.Compound;
import com.insilico.application.insilicotools.data.Project;
import com.insilico.application.insilicotools.data.Status;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import openeye.oechem.OEGraphMol;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

public class MyCompoundModifyDialog extends StandardDialog {
    Compound myCompound;
    JTextField nameField = new JTextField(10);
    JComboBox<Status> statusCB = new JComboBox<>(FrontierDAO.getInstance().getMy_status());
    MSketchPane sketcher = MarvinFactory.createCompoundSketcher();
    JButton okBtn;
    JButton cancelBtn;
    Project project;

    public Compound getMyCompound() {
        return myCompound;
    }

    public MyCompoundModifyDialog(Project project1) {
        this.project = project1;
        okBtn = new JButton("OK");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(nameField.getText().isEmpty()){
                    JOptionPane.showMessageDialog(MyCompoundModifyDialog.this,"No name specified.");
                    return;
                }
                if(sketcher.getMol().isEmpty()){
                    JOptionPane.showMessageDialog(MyCompoundModifyDialog.this,"No structure specified.");
                    return;
                }
                Molecule mol = sketcher.getMol();
                String mol_name =nameField.getText();
                OEGraphMol oeGraphMol = OEChemFunc.getInstance().convertChemAxonMol(mol);
                String molStr = ChemFunc.getMolString(oeGraphMol);
                byte[] cdx = new byte[0];
                try {
                    cdx = MolExporter.exportToBinFormat(mol, "cdx");
                } catch (IOException e1) {
                    e1.printStackTrace();
                }

                myCompound = new Compound(-1,nameField.getText(),sketcher.getMol("sdf"),sketcher.getMol("smiles:u,a"),cdx, InSlilicoPanel.getInstance().getUserName(), new java.sql.Date(System.currentTimeMillis()), project.getProject_id(),statusCB.getItemAt(statusCB.getSelectedIndex()).getStatus_id(),0,0);                myCompound.setStatus_id(statusCB.getItemAt(statusCB.getSelectedIndex()).getStatus_id());
                myCompound.setMol(molStr);
                myCompound.setName(mol_name);
                myCompound.setChanged(true);
                MyCompoundModifyDialog.this.dispose();
            }
        });
        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MyCompoundModifyDialog.this.dispose();
            }
        });
        setSize(new Dimension(1024,768));
    }

    public MyCompoundModifyDialog(Compound compound, Project project) {
        this.project = project;
        myCompound = compound;
        if(myCompound !=null){
            Molecule mol = OEChemFunc.getInstance().convertOEChemMol(myCompound.getPropertyMol().getMol());
            sketcher.setMol(mol);
            nameField.setText(myCompound.getName());
            statusCB.setSelectedItem(myCompound.getStatus());
        }
        okBtn = new JButton("OK");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(nameField.getText().isEmpty()){
                    JOptionPane.showMessageDialog(MyCompoundModifyDialog.this,"No name specified.");
                    return;
                }
                if(sketcher.getMol().isEmpty()){
                    JOptionPane.showMessageDialog(MyCompoundModifyDialog.this,"No structure specified.");
                    return;
                }
                myCompound.setStatus_id(statusCB.getItemAt(statusCB.getSelectedIndex()).getStatus_id());
                myCompound.setMol(sketcher.getMol("sdf"));
                myCompound.setName(nameField.getText());
                MyCompoundModifyDialog.this.dispose();
            }
        });
        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MyCompoundModifyDialog.this.dispose();
            }
        });
        setSize(new Dimension(1024,768));
    }

    @Override
    public JComponent createBannerPanel() {
        return null;
    }

    @Override
    public JComponent createContentPanel() {
        JPanel p = new JPanel(new BorderLayout());
        p.add(sketcher, BorderLayout.CENTER);
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Name:"));
        p1.add(nameField);
        p1.add(new JLabel("Status"));
        p1.add(statusCB);
        p.add(p1,BorderLayout.SOUTH);
        return p;
    }

    @Override
    public ButtonPanel createButtonPanel() {
        ButtonPanel p = new ButtonPanel();
        p.addButton(okBtn,ButtonPanel.AFFIRMATIVE_BUTTON);
        p.addButton(cancelBtn,ButtonPanel.CANCEL_BUTTON);
        return p;
    }

    public static void main(String[] args) {

    }
}
