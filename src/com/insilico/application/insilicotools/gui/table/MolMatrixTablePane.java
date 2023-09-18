package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.util.StatFunc;
import openeye.oechem.OEFormat;
import openeye.oechem.oechem;
import openeye.oechem.oemolostream;
import org.jdesktop.swingx.JXTable;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.*;

/**
 * Created by jfeng1 on 2/17/16.
 */
public class MolMatrixTablePane extends JPanel {
    MatrixCellTable cellTable;
    MatrixMolTableModel molTableModel;
    JComboBox columnCountCB;
    JSlider columnWidthSlider;
    JButton saveAsSdfBtn;
    File currentDirectory;
    DesignProgressMonitor progressMonitor = new DesignProgressMonitor(this,"Progress","Progress",0,100);
    boolean simplified = false;
    TitledBorder border;
    JPanel btnPanel;

    public MolMatrixTablePane(MatrixMolTableModel molTableModel) {
        this(molTableModel,false);
    }

    public MolMatrixTablePane(final MatrixMolTableModel molTableModel, boolean simplified) {
        super(new BorderLayout());
        border = new TitledBorder("Molecules (0)");
        this.molTableModel = molTableModel;
        this.simplified = simplified;
        this.molTableModel.addTableModelListener(new TableModelListener() {
            @Override
            public void tableChanged(TableModelEvent e) {
                border.setTitle(String.format("Molecules(%d)",molTableModel.getPropertyMolecules().size()));
                fixMatrixTableFormat(cellTable,molTableModel,columnWidthSlider.getValue());
                MolMatrixTablePane.this.repaint();
            }
        });
        cellTable = new MatrixCellTable(molTableModel);
        add(new JScrollPane(cellTable),BorderLayout.CENTER);
        btnPanel = buildBtnPanel();
        add(btnPanel, BorderLayout.SOUTH);
        if(!simplified) {
            add(buildTopPanel(), BorderLayout.NORTH);
        }
        setBorder(border);
    }

    public JPanel getBtnPanel() {
        return btnPanel;
    }

    private void fixMatrixTableFormat(JXTable matrixTable, MatrixMolTableModel matrixMolTableModel, int columnWidth){
        matrixTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        for(int i=0;i<matrixMolTableModel.getColumnCount();i++){
            matrixTable.getColumnModel().getColumn(i).setPreferredWidth(150+columnWidth);
        }
        matrixTable.setRowHeight(150+columnWidth);
    }

    private JPanel buildTopPanel(){
        JPanel topPanel = new JPanel();
        JButton selectAllBtn = new JButton("Select All");
        selectAllBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                cellTable.selectAll();
            }
        });
        JButton invertSelectionBtn = new JButton("Invert Selection");
        invertSelectionBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                cellTable.invertSelection();
            }
        });
        JButton unselectAllBtn = new JButton("Unselect All");
        unselectAllBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                cellTable.unselectAll();
            }
        });
        JButton deleteSelectBtn = new JButton("Delete selected");
        deleteSelectBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int option = JOptionPane.showConfirmDialog(MolMatrixTablePane.this, "Selected molecules will be deleted, are you sure ?");
                if(option!=JOptionPane.YES_OPTION){
                    return;
                }
                Vector<PropertyMolecule> molecules = molTableModel.getPropertyMolecules();
                Vector<PropertyMolecule> molToBeRemoved = new Vector<PropertyMolecule>();
                for(PropertyMolecule mol:molecules){
                    if(mol.isSelected()){
                        molToBeRemoved.add(mol);
                    }
                }
                for(PropertyMolecule mol:molToBeRemoved){
                    molecules.remove(mol);
                }
                molTableModel.fireTableDataChanged();
            }
        });

        final JTextField numCmpdsField = new JTextField(10);
        JButton diverseBtn = new JButton("Select ");
        diverseBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    int numCmpds = Integer.parseInt(numCmpdsField.getText());
                    if(numCmpds>molTableModel.getPropertyMolecules().size()){
                        throw new NumberFormatException("Too many molecules.");
                    }
                    Vector<PropertyMolecule> diverseMolecules = StatFunc.getDiverseMolecules(molTableModel.getPropertyMolecules(), numCmpds);
                    for(PropertyMolecule mol:diverseMolecules){
                        mol.setIsSelected(true);
                    }
                    molTableModel.fireTableDataChanged();
                } catch (NumberFormatException e1) {
                    JOptionPane.showMessageDialog(MolMatrixTablePane.this,e1.getMessage());
                }
            }
        });

        JButton filterBtn = new JButton("Filter");
        filterBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Vector<PropertyMolecule> molecules = molTableModel.getPropertyMolecules();
                Vector<PropertyMolecule> molToBeRemoved = new Vector<PropertyMolecule>();
                for(PropertyMolecule mol:molecules){
                    if(mol.getNumAromaticRings()==0){
                        molToBeRemoved.add(mol);
                    }
                }
                for(PropertyMolecule mol:molToBeRemoved){
                    molecules.remove(mol);
                }
                molTableModel.fireTableDataChanged();
            }
        });

        topPanel.add(selectAllBtn);
        topPanel.add(invertSelectionBtn);
        topPanel.add(unselectAllBtn);
        topPanel.add(deleteSelectBtn);
//        if (!simplified) {
//            topPanel.add(diverseBtn);
//            topPanel.add(numCmpdsField);
//            topPanel.add(new JLabel(" diverse compounds."));
//        }
        return topPanel;
    }

    private JPanel buildBtnPanel(){
        JPanel btnPanel = new JPanel();
        columnCountCB = new JComboBox(new Integer[]{2,3,4,5,6,7,8});
        columnCountCB.setSelectedItem(new Integer(4));
        columnCountCB.setMaximumSize(new Dimension(100,20));
        columnCountCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                molTableModel.setColumnCount((Integer)columnCountCB.getSelectedItem());
                fixMatrixTableFormat(cellTable,molTableModel,columnWidthSlider.getValue());
            }
        });

        columnWidthSlider = new JSlider(0,250,0);
        columnWidthSlider.setMaximumSize(new Dimension(200,20));
        columnWidthSlider.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                fixMatrixTableFormat(cellTable,molTableModel,columnWidthSlider.getValue());
            }
        });
        saveAsSdfBtn = new JButton("Save As Sdf");
        saveAsSdfBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                final Vector<PropertyMolecule> propertyMolecules = molTableModel.getPropertyMolecules();
                if(propertyMolecules.size()==0){
                    JOptionPane.showMessageDialog(MolMatrixTablePane.this,"No molecules available.");
                    return;
                }
                FileUtil.saveToFile(currentDirectory, new FileNameExtensionFilter("SDF file","sdf"), new FileFunctor() {
                    @Override
                    public void execute(final File file) {
                        SwingWorker sw = new SwingWorker() {
                            @Override
                            protected Object doInBackground() throws Exception {
                                oemolostream ofs = new oemolostream();
                                ofs.SetFormat(OEFormat.SDF);
                                ofs.open(file.getAbsolutePath());
                                int progress = 0;
                                for(PropertyMolecule mol:propertyMolecules){
                                    progress ++;
                                    Vector v = new Vector();
                                    v.add(String.format("Saving molecule No. %d",progress));
                                    v.add(100*progress/propertyMolecules.size());
                                    publish(v);
                                    oechem.OEWriteMolecule(ofs,mol.getMol());
                                }
                                ofs.close();
                                return null;
                            }

                            @Override
                            protected void process(java.util.List chunks) {
                                Vector v = (Vector) chunks.get(chunks.size()-1);
                                String note = (String)v.get(0);
                                int progress = (Integer) v.get(1);
                                progressMonitor.setNote(note);
                                progressMonitor.setProgress(progress);
                            }

                            @Override
                            protected void done() {
                                progressMonitor.close();
                                try {
                                    get();
                                    JOptionPane.showMessageDialog(MolMatrixTablePane.this, "SDF file saved.");
                                } catch (Exception e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(MolMatrixTablePane.this,e1.getMessage());
                                }finally {

                                }
                            }
                        };
                        sw.execute();

                    }
                });
            }
        });
        btnPanel.add(new JLabel("Column Count:"));
        btnPanel.add(columnCountCB);
        btnPanel.add(new JSeparator());
        btnPanel.add(new JLabel("Column Width:"));
        btnPanel.add(columnWidthSlider);
        btnPanel.add(saveAsSdfBtn);
        return btnPanel;
    }
}
