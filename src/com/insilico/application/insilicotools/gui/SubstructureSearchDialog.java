package com.insilico.application.insilicotools.gui;

import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.insilico.application.insilicotools.gui.filter.smarts.SmartsPattern;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.gui.widget.FloatTextField;
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
public class SubstructureSearchDialog extends StandardDialog{
    MSketchPane sketcher;
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;
    private int maxNumHits;
    JTextField numHitsField;
    boolean useNumHits = true;
    boolean useProperty = true;
    FloatTextField mwField;
    FloatTextField clogpField;
    FloatTextField psaField;


    public SubstructureSearchDialog(Frame frame, String title, boolean useNumHits, boolean useProperty) {
        super(frame,title);
        setModal(true);
        setSize(new Dimension(800,800));
        this.useNumHits = useNumHits;
        this.useProperty = useProperty;
    }

    public SubstructureSearchDialog(Frame frame, String title) {
        this(frame,title,true,true);
    }

    public SubstructureSearchDialog(JDialog dialog, String title){
        this(dialog,title,true,true);
    }
    public SubstructureSearchDialog(JDialog dialog, String title, boolean useNumHits, boolean useProperty) {
        super(dialog,title);
        setModal(true);
        setSize(new Dimension(800,800));
        this.useNumHits = useNumHits;
        this.useProperty = useProperty;
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public String getSmarts(){
        if(!sketcher.getMol().isEmpty()){
            return sketcher.getMol("smarts:ah");
        }
        return null;
    }

    public float getCLogPLimit(){
        if(clogpField!=null&&!clogpField.getText().trim().isEmpty()){
            float v = Float.parseFloat(clogpField.getText());
            if(v>0){
                return v;
            }
        }
        return -1;
    }

    public float getMWLimit(){
        if(mwField!=null&&!mwField.getText().trim().isEmpty()){
            float v = Float.parseFloat(mwField.getText());
            if(v>0){
                return v;
            }
        }
        return -1;
    }

    public float getPSALimit(){
        if(psaField!=null&&!psaField.getText().trim().isEmpty()){
            float v = Float.parseFloat(psaField.getText());
            if(v>0){
                return v;
            }
        }
        return -1;
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
                    JOptionPane.showMessageDialog(SubstructureSearchDialog.this,"No CY Number Available.");
                    return;
                }
                String formattedCYNumber = ChemFunc.formatCYNumber(bioField.getText());
                if(formattedCYNumber==null){
                    JOptionPane.showMessageDialog(SubstructureSearchDialog.this,"Illegal CY Number Available.");
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
        JPanel bottomPanel = new JPanel(new GridLayout(2,1));
        if(useNumHits) {
            JPanel p1 = new JPanel();
            final JComboBox smartsCB = new JComboBox(SmartsPattern.getPredefinedSmartsPattern());
            smartsCB.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    SmartsPattern smartsPattern = (SmartsPattern) smartsCB.getSelectedItem();
                    sketcher.setMol(smartsPattern.getSmarts());
                    okBtn.setEnabled(true);
                }
            });
            p1.add(new JPopupMenu.Separator());
            p1.add(new JLabel("Predefined Smarts:"));
            p1.add(smartsCB);

            numHitsField = new JTextField(5);
            //numHitsField.setText(""+DEFAULT_NUM_HITS);
            p1.add(new JLabel("Max Hits:"));
            p1.add(numHitsField);
            bottomPanel.add(p1);
        }

        if(useProperty) {
            JPanel p2 = new JPanel();
            clogpField = new FloatTextField(5);
            mwField = new FloatTextField(5);
            psaField = new FloatTextField(5);
            p2.add(new JLabel("CLogP < "));
            p2.add(clogpField);
            p2.add(new JLabel("MW < "));
            p2.add(mwField);
            p2.add(new JLabel("PSA < "));
            p2.add(psaField);
            bottomPanel.add(p2);
        }
        if(useProperty||useNumHits) {
            mainPanel.add(bottomPanel, BorderLayout.SOUTH);
        }
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
                SubstructureSearchDialog.this.setVisible(false);
            }
        });

        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                SubstructureSearchDialog.this.setVisible(false);
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
        SubstructureSearchDialog dialog = new SubstructureSearchDialog(new JFrame(),"");
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setVisible(true);
    }

    public int getMaxNumHits() {
        if(useNumHits) {
            try {
                maxNumHits = Integer.parseInt(numHitsField.getText());
                return maxNumHits;
            } catch (NumberFormatException e) {
                return -1;
            }
        }else{
            return -1;
        }
    }
}
