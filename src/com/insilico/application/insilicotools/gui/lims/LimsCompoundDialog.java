package com.insilico.application.insilicotools.gui.lims;

import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.database.LimsDAO;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.gui.widget.FloatTextField;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import openeye.oechem.OEGraphMol;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import java.sql.SQLException;

public class LimsCompoundDialog extends StandardDialog {
    JTextField nameField = new JTextField(10);
    FloatTextField theoreticalMassField = new FloatTextField(5);
    FloatTextField amountSubmittedField = new FloatTextField(5);
    JComboBox compatibilityWithTFACB = new JComboBox(new String[]{"yes","no"});
    JComboBox<String> unitCB = new JComboBox<>(new String[]{"mg","g","ug"});
    JComboBox<String> projectCB = new JComboBox<>(LimsDAO.getInstance().getMyProjectNames());
    MSketchPane sketcher = MarvinFactory.createCompoundSketcher();
    JButton okBtn;
    JButton cancelBtn;
    LimsMolecule targetMol;
    boolean isCommitted = false;


    public boolean isCommitted() {
        return isCommitted;
    }

    public LimsCompoundDialog() {

        okBtn = new JButton("OK");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(nameField.getText().isEmpty()){
                    JOptionPane.showMessageDialog(LimsCompoundDialog.this,"No name specified.");
                    return;
                }
                if(sketcher.getMol().isEmpty()){
                    JOptionPane.showMessageDialog(LimsCompoundDialog.this,"No structure specified.");
                    return;
                }
                String project = (String)projectCB.getSelectedItem();
                Molecule mol = sketcher.getMol();
                String compound_id =nameField.getText();
                OEGraphMol oeGraphMol = OEChemFunc.getInstance().convertChemAxonMol(mol);
                String compatibility_with_tfa = (String)compatibilityWithTFACB.getSelectedItem();
                String amount_submitted = amountSubmittedField.getText();
                if(amount_submitted.isEmpty()){
                    amount_submitted = null;
                }else{
                    amount_submitted=amount_submitted+" mg";
                }
                String theory_mass = theoreticalMassField.getText() + " "+ unitCB.getSelectedItem();
                targetMol = new LimsMolecule(new PropertyMolecule(oeGraphMol), -1, compound_id, theory_mass, project, InSlilicoPanel.getInstance().getUserName(), compatibility_with_tfa, amount_submitted, System.currentTimeMillis());
                try {
                    LimsDAO.getInstance().insertLimsCompound(targetMol);
                    isCommitted = true;
                    LimsCompoundDialog.this.dispose();
                } catch (SQLException e1) {
                    e1.printStackTrace();
                }
            }
        });
        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                LimsCompoundDialog.this.dispose();
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
        JPanel btmPanel = new JPanel(new GridLayout(2,1));
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Compound ID:"));
        p1.add(nameField);
        p1.add(new JLabel("Theoretical Mass:"));
        p1.add(theoreticalMassField);
        p1.add(unitCB);
        p1.add(new JLabel("Compatibility With TFA:"));
        p1.add(compatibilityWithTFACB);
        p1.add(new JLabel("Amount Submitted (mg):"));
        p1.add(amountSubmittedField);
        JPanel p11 = new JPanel();
        p11.add(new JLabel("Project:"));
        p11.add(projectCB);
        btmPanel.add(p1);
        btmPanel.add(p11);

        p.add(btmPanel,BorderLayout.SOUTH);
        JPanel p2 = new JPanel();
        p2.add(new JLabel("BIO-"));
        final JTextField biofield = new JTextField(10);
        p2.add(biofield);
        JButton loadBtn = new JButton("Load");
        loadBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String cynumber = ChemFunc.formatCYNumber(biofield.getText());
                if(cynumber!=null){
                    try {
                        OEGraphMol cy_mol = InHouseCollectionDAO.getInstance().getMolFromCYNumber(cynumber);
                        if(cy_mol!=null) {
                            sketcher.setMol(OEChemFunc.getInstance().convertOEChemMol(cy_mol));
                        }
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(null,e1.getMessage());
                    }
                }else{
                    JOptionPane.showMessageDialog(LimsCompoundDialog.this,"Wrong BioNumber Provided.");
                }

            }
        });
        p2.add(loadBtn);
        p.add(p2,BorderLayout.NORTH);
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
