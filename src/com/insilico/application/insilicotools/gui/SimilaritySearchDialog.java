package com.insilico.application.insilicotools.gui;

import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
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
public class SimilaritySearchDialog extends StandardDialog{
    MSketchPane sketcher;
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;
    private int maxNumHits;
    private float simCutoff;
    JTextField numHitsField;
    JTextField simCutoffField;

    public SimilaritySearchDialog(Frame frame, String title) {
        super(frame,title);
        setModal(true);
        setSize(new Dimension(800,800));
    }

    public SimilaritySearchDialog(JDialog dialog, String title) {
        super(dialog,title);
        setModal(true);
        setSize(new Dimension(800,800));
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public Molecule getMolecule(){
        return sketcher.getMol();
    }

    @Override
    public JComponent createBannerPanel() {
        JPanel p = new JPanel();
        final JTextField bioField = new JTextField(10);
        final JButton btn = new JButton("Load");
        btn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(bioField.getText().isEmpty()){
                    JOptionPane.showMessageDialog(SimilaritySearchDialog.this,"No CY Number Available.");
                    return;
                }
                String formattedCYNumber = ChemFunc.formatCYNumber(bioField.getText());
                if(formattedCYNumber==null){
                    JOptionPane.showMessageDialog(SimilaritySearchDialog.this,"Illegal CY Number Available.");
                    return;
                }
                try {
                    OEGraphMol molFromCYNumber = InHouseCollectionDAO.getInstance().getMolFromCYNumber(formattedCYNumber);
                    if(molFromCYNumber!=null) {
                        bioField.setText(formattedCYNumber);
                        Molecule mol = OEChemFunc.getInstance().convertOEChemMol(molFromCYNumber);
                        mol.setName(formattedCYNumber);
                        sketcher.setMol(mol);
                        okBtn.setEnabled(true);
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
        p.add(new JLabel("CY-NUMBER:"));
        p.add(bioField);
        p.add(btn);

        return p;
    }

    @Override
    public JComponent createContentPanel() {
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBorder(BorderFactory.createTitledBorder("Sketch a molecule:"));
        sketcher = MarvinFactory.createCompoundSketcher();
        sketcher.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(evt.getPropertyName().equals("mol")){
                    if(sketcher.getMol().isEmpty()){
                        okBtn.setEnabled(false);
                    }else{
                        okBtn.setEnabled(true);
                    }
                }
            }
        });
        mainPanel.add(sketcher,BorderLayout.CENTER);
        JPanel p = new JPanel();

        simCutoffField = new JTextField(5);
        simCutoffField.setText("0.7");
        p.add(new JLabel("Sim. Cutoff:"));
        p.add(simCutoffField);

        numHitsField = new JTextField(5);
        //numHitsField.setText(""+DEFAULT_NUM_HITS);
        p.add(new JLabel("Max Hits:"));
        p.add(numHitsField);
        mainPanel.add(p,BorderLayout.SOUTH);
//        mainPanel.add(bottomPanel,BorderLayout.SOUTH);
        return mainPanel;
    }

    @Override
    public ButtonPanel createButtonPanel() {
        ButtonPanel buttonPanel = new ButtonPanel(SwingConstants.CENTER);

        okBtn = new JButton("OK");
        okBtn.setEnabled(false);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(!sketcher.getMol().isEmpty()){
                    isCommitted = true;
                }else{
                    isCommitted = false; //extra layer of protection from empty molecule
                }
                SimilaritySearchDialog.this.setVisible(false);
            }
        });

        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                SimilaritySearchDialog.this.setVisible(false);
            }
        });
        buttonPanel.addButton(okBtn,ButtonPanel.AFFIRMATIVE_BUTTON);
        buttonPanel.addButton(cancelBtn,ButtonPanel.CANCEL_BUTTON);
        return buttonPanel;
    }

    public void setVisible(boolean b) {
        super.setVisible(b);
        if (b) {
            if(sketcher!=null) {
                sketcher.requestFocus();
            }
        }
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        SimilaritySearchDialog dialog = new SimilaritySearchDialog(f,"");
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setVisible(true);
    }

    public float getSimCutoff(){
        try{
            simCutoff = Float.parseFloat(simCutoffField.getText());
            return simCutoff;
        }catch (NumberFormatException e){
            return 0.7F;
        }

    }

    public int getMaxNumHits() {
        try {
            maxNumHits = Integer.parseInt(numHitsField.getText());
            return maxNumHits;
        } catch (NumberFormatException e) {
            return -1;
        }
    }
}
