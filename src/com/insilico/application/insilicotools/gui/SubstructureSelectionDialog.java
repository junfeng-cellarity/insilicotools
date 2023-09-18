package com.insilico.application.insilicotools.gui;

import chemaxon.calculations.clean.Cleaner;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.marvin.paint.DispOptConsts;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import openeye.oechem.OEGraphMol;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.util.Vector;

/**
 * Created by jfeng1 on 1/4/16.
 */
public class SubstructureSelectionDialog extends StandardDialog{
    MSketchPane sketcher;
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;
    Molecule templateMol;
    static boolean updating = false;
    String smarts;
    String atomIdsStr;

    public SubstructureSelectionDialog() {
        super();
        init();

    }

    private void init() {
        setModal(true);
        setSize(new Dimension(800,800));
        sketcher = MarvinFactory.createCompoundSketcher();
        sketcher.setImplicitH(DispOptConsts.IMPLICITH_OFF_S);
        sketcher.setSketchMode(MSketchPane.SM_SELECT_LASSO);
        okBtn = new JButton("OK");
        cancelBtn = new JButton("Cancel");
    }

    public SubstructureSelectionDialog(JFrame parentFrame) {
        super(parentFrame);
        init();
    }


    public boolean isCommitted() {
        return isCommitted;
    }

    public Molecule getMolecule(){
        Molecule mol = sketcher.getMol();
        return mol;
    }

    public String getSmarts() {
        return smarts;
    }


    public void setMolecule(Molecule mol){
        if(mol!=null) {
            sketcher.setMol(mol);
            sketcher.firePropertyChange("selectionChanged",true,false);
        }
    }

    public void setPropertyMolecule(PropertyMolecule mol){
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
        }
    }

    @Override
    public JComponent createBannerPanel() {
        JPanel p = new JPanel();
        p.add(new JLabel("Select a substructure as constraint. Must less than 30 heavy atoms."));
        return p;
    }

    @Override
    public JComponent createContentPanel() {
        JPanel mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBorder(BorderFactory.createTitledBorder("Select substructure:"));
        sketcher.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(evt.getPropertyName().equals("mol")){
                    if(updating){
                        return;
                    }
                    updating = true;
                    if(templateMol!=null) {
                        sketcher.setMol(templateMol);
                    }
                    if(sketcher.getMol().isEmpty()){
                        okBtn.setEnabled(false);
                    }else{
                        okBtn.setEnabled(true);
                    }
                    updating = false;
                }else if(evt.getPropertyName().equals("selectionChanged")){

                    extractSelectedAtoms();

                }
            }
        });
        mainPanel.add(sketcher,BorderLayout.CENTER);
//        mainPanel.add(bottomPanel,BorderLayout.SOUTH);
        return mainPanel;
    }

    private void extractSelectedAtoms() {
        Molecule m = getMolecule().clone();
        Vector<MolAtom> selectedAtoms = new Vector<MolAtom>();
        StringBuilder sb = new StringBuilder();
        for(MolAtom a:m.getAtomArray()){
            if(a.isSelected()){
                if(!selectedAtoms.isEmpty()){
                    sb.append(", ");
                }
                selectedAtoms.add(a);
                sb.append(m.indexOf(a)+1);
            }
        }
        for(MolBond b:m.getBondArray()){
            if(selectedAtoms.contains(b.getAtom1())&&selectedAtoms.contains(b.getAtom2())){
                continue;
            }
            m.removeBond(b);
        }
        for(MolAtom a:m.getAtomArray()){
            if(!selectedAtoms.contains(a)){
                m.removeAtom(a);
            }
        }
        try {
            int hvyAtomCount = 0;
            for(MolAtom a:m.getAtomArray()){
                if(a.getAtno()>1){
                    hvyAtomCount += 1;
                }
            }
            if(hvyAtomCount>29){
                throw new IOException("Substructure must have less than 30 heavy atoms.");
            }
            if(hvyAtomCount<5){
                throw new IOException("Substructure must have more than 4 heavy atoms.");
            }

            smarts = MolExporter.exportToFormat(m,"smiles:a");
            if(smarts.contains(".")){
                smarts = null;
                throw new IOException("Must be connected substructure.");
            }
            atomIdsStr = sb.toString();
            System.out.println(atomIdsStr);
            System.out.println(smarts);
            okBtn.setEnabled(true);
        } catch (IOException e) {
            okBtn.setEnabled(false);
            System.out.println(e.getMessage());
        }
    }

    public String getAtomIdsStr() {
        return atomIdsStr;
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
                SubstructureSelectionDialog.this.setVisible(false);
            }
        });

        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                SubstructureSelectionDialog.this.setVisible(false);
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
        SubstructureSelectionDialog dialog = new SubstructureSelectionDialog();
        String s = "BIO-0704342\n" +
                "  -OEChem-10041610313D\n" +
                "\n" +
                " 38 40  0     0  0  0  0  0  0999 V2000\n" +
                "    9.7058    2.1007   -4.2335 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    8.9385    3.1156   -3.6743 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.1120    0.7948   -0.3873 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.0496   -1.3108   -0.2349 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0387    1.4326    0.2455 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.1241   -0.6764    0.3979 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    9.2130    0.7990   -4.2468 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.4251   -0.8216   -1.7975 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    7.7012    2.7874   -3.1471 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.0847   -0.5866   -0.6374 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    4.0940   -1.9660   -2.3340 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.0746    0.6935    0.6352 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.2091   -1.2463   -1.2988 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    7.9579    0.5509   -3.6959 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.2221   -3.0208   -2.1192 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.3957   -2.0843   -2.9770 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.2835   -4.4675   -2.4237 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    7.4027   -0.8577   -3.6998 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    7.1924    1.5353   -3.1433 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    6.0958   -0.9132   -3.0875 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.8391   -3.1513   -3.3916 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.0905   -2.5906   -1.4980 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -2.3933    1.4775    1.4128 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   10.6819    2.3214   -4.6567 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    9.2966    4.1387   -3.6499 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.9689    1.3959   -0.6811 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.1075   -2.3835   -0.4116 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0807    2.5029    0.4294 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.9902   -1.2591    0.6996 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    9.8007   -0.0043   -4.6801 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.7871    0.1958   -1.7765 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    7.0671    3.5477   -2.6999 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.3773   -4.9827   -2.0895 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.3841   -4.6284   -3.5014 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    4.1383   -4.9292   -1.9203 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    7.3141   -1.2222   -4.7291 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    8.0701   -1.5246   -3.1429 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    5.7373   -0.0300   -2.7482 H   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  1  0  0  0  0\n" +
                "  1  7  2  0  0  0  0\n" +
                "  2  9  2  0  0  0  0\n" +
                "  3  5  1  0  0  0  0\n" +
                "  3 10  2  0  0  0  0\n" +
                "  4  6  2  0  0  0  0\n" +
                "  4 10  1  0  0  0  0\n" +
                "  5 12  2  0  0  0  0\n" +
                "  6 12  1  0  0  0  0\n" +
                "  7 14  1  0  0  0  0\n" +
                "  8 11  1  0  0  0  0\n" +
                "  8 13  2  0  0  0  0\n" +
                "  9 19  1  0  0  0  0\n" +
                " 10 13  1  0  0  0  0\n" +
                " 11 15  2  0  0  0  0\n" +
                " 11 16  1  0  0  0  0\n" +
                " 12 23  1  0  0  0  0\n" +
                " 13 22  1  0  0  0  0\n" +
                " 14 18  1  0  0  0  0\n" +
                " 14 19  2  0  0  0  0\n" +
                " 15 17  1  0  0  0  0\n" +
                " 15 22  1  0  0  0  0\n" +
                " 16 20  1  0  0  0  0\n" +
                " 16 21  2  0  0  0  0\n" +
                " 18 20  1  0  0  0  0\n" +
                "  1 24  1  0  0  0  0\n" +
                "  2 25  1  0  0  0  0\n" +
                "  3 26  1  0  0  0  0\n" +
                "  4 27  1  0  0  0  0\n" +
                "  5 28  1  0  0  0  0\n" +
                "  6 29  1  0  0  0  0\n" +
                "  7 30  1  0  0  0  0\n" +
                "  8 31  1  0  0  0  0\n" +
                "  9 32  1  0  0  0  0\n" +
                " 17 33  1  0  0  0  0\n" +
                " 17 34  1  0  0  0  0\n" +
                " 17 35  1  0  0  0  0\n" +
                " 18 36  1  0  0  0  0\n" +
                " 18 37  1  0  0  0  0\n" +
                " 20 38  1  0  0  0  0\n" +
                "M  END\n" +
                "> <Energy>\n" +
                "13.1133\n" +
                "\n" +
                "$$$$\n";
        Molecule m = null;
        try {
            m = MolImporter.importMol(s);
        } catch (MolFormatException e) {
            e.printStackTrace();
            return;
        }
        OEGraphMol mol = OEChemFunc.getInstance().convertChemAxonMol(m);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setPropertyMolecule(new PropertyMolecule(mol));
        dialog.setVisible(true);

    }
}
