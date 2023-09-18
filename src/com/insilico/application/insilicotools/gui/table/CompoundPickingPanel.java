package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.filter.FilterDialog;
import com.insilico.application.insilicotools.gui.util.ImageUtil;
import com.jidesoft.swing.JideSplitPane;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

/**
 * Created by jfeng1 on 3/30/16.
 */
public class CompoundPickingPanel extends JPanel {
    MolMatrixTablePane leftTable;
    MolMatrixTablePane rightTable;
    MatrixMolTableModel leftTableModel;
    MatrixMolTableModel rightTableModel;
    Vector<PropertyMolecule> unselectedMolecules = new Vector<PropertyMolecule>();
    Vector<PropertyMolecule> selectedMolecules = new Vector<PropertyMolecule>();

    public MatrixMolTableModel getRightTableModel() {
        return rightTableModel;
    }

    public MatrixMolTableModel getLeftTableModel() {
        return leftTableModel;
    }

    private void initializeComponents() {
        JPanel buttonPanel = buildButtonPanel();
        JideSplitPane splitPane = new JideSplitPane(JideSplitPane.HORIZONTAL_SPLIT);
        leftTableModel = new MatrixMolTableModel(unselectedMolecules);
        leftTable = new MolMatrixTablePane(leftTableModel);
        rightTableModel = new MatrixMolTableModel(selectedMolecules);
        rightTable = new MolMatrixTablePane(rightTableModel);
        splitPane.add(leftTable);
        splitPane.add(rightTable);
        add(splitPane, BorderLayout.CENTER);
        add(buttonPanel,BorderLayout.SOUTH);
    }

    public CompoundPickingPanel(Vector<PropertyMolecule> molecules1, Vector<PropertyMolecule> molecules2){
        super(new BorderLayout());
        initializeComponents();
        if(molecules1==null){
            unselectedMolecules = new Vector<PropertyMolecule>();
        }else{
            unselectedMolecules = molecules1;
        }
        if(molecules2==null) {
            selectedMolecules = new Vector<PropertyMolecule>();
        }else{
            selectedMolecules = molecules2;
        }
        initializeComponents();
    }

    private JPanel buildButtonPanel(){
        JPanel p = new JPanel();
        p.setPreferredSize(new Dimension(0,100));
        ImageIcon filterIcon = ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("filter-icon.png")));
        final JButton filterBtn = new JButton("Filter",filterIcon);
        filterBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                PropertyMolecule[] initialMolecules = unselectedMolecules.toArray(new PropertyMolecule[unselectedMolecules.size()]);
                Container topLevelAncestor = CompoundPickingPanel.this.getTopLevelAncestor();
                FilterDialog filterDialog = null;
                if(topLevelAncestor instanceof JFrame){
                    filterDialog = new FilterDialog((JFrame) topLevelAncestor,initialMolecules);
                }else{
                    filterDialog = new FilterDialog((JDialog) topLevelAncestor,initialMolecules);
                }

                filterDialog.setLocationRelativeTo(CompoundPickingPanel.this);
                filterDialog.setVisible(true);
                if(filterDialog.isCommitted()){
                    PropertyMolecule[] molInShoppingCart = filterDialog.getSelectedMolecules();
                    if(molInShoppingCart!=null&&molInShoppingCart.length>0) {
                        unselectedMolecules.clear();
                        for (PropertyMolecule mol : molInShoppingCart) {
                            unselectedMolecules.add(mol);
                        }
                        leftTableModel.fireTableDataChanged();
                    }
                }
            }
        });
        ImageIcon insideIcon = ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("Arrow-Inside-icon.png")));
        JButton addToSelectBtn = new JButton("Add Selected",insideIcon);
        addToSelectBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Vector<PropertyMolecule> v = new Vector<PropertyMolecule>();
                for(PropertyMolecule m:unselectedMolecules){
                    if(m.isSelected()){
                        v.add(m);
                    }
                }
                for(PropertyMolecule m:v){
                    if(!selectedMolecules.contains(m)) {
                        m.setIsSelected(false);
                        selectedMolecules.add(m);
                    }
                    unselectedMolecules.remove(m);
                }
                leftTableModel.fireTableDataChanged();
                rightTableModel.fireTableDataChanged();
            }
        });
        ImageIcon outsiteIcon = ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("Arrow-Outside-icon.png")));
        JButton removeToSelectBtn = new JButton("Remove Selected",outsiteIcon);
        removeToSelectBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Vector<PropertyMolecule> v = new Vector<PropertyMolecule>();
                for(PropertyMolecule m:selectedMolecules){
                    if(m.isSelected()){
                        v.add(m);
                    }
                }
                for(PropertyMolecule m:v){
                    if(!unselectedMolecules.contains(m)) {
                        m.setIsSelected(false);
                        unselectedMolecules.add(m);
                    }
                    selectedMolecules.remove(m);
                }
                leftTableModel.fireTableDataChanged();
                rightTableModel.fireTableDataChanged();
            }
        });
        p.add(addToSelectBtn);
        p.add(filterBtn);
        p.add(removeToSelectBtn);
        return p;
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        Vector<PropertyMolecule> mols = new Vector<PropertyMolecule>();
        oemolistream ifs = new oemolistream();
        ifs.open("/Users/jfeng1/00Demo/demo.sdf");
        OEGraphMol mol = new OEGraphMol();
        while(oechem.OEReadMolecule(ifs,mol)){
            mols.add(new PropertyMolecule(mol));
        }

        CompoundPickingPanel p = new CompoundPickingPanel(mols, new Vector<PropertyMolecule>());
        f.getContentPane().add(p);
        f.setSize(new Dimension(1024,768));
        f.setVisible(true);
    }


}
