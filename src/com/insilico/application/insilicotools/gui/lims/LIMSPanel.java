package com.insilico.application.insilicotools.gui.lims;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.LimsDAO;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.table.MolTable2D;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.jidesoft.swing.JideSplitPane;
import openeye.oechem.*;
import org.jdesktop.swingx.JXTable;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.TableCellEditor;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.sql.SQLException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 3/23/17.
 */
public class LIMSPanel extends JPanel {
    BatchTableModel batchTableModel;
    MolTable2D molTable2D;
    Vector<LimsMolecule> limsMolecules;
    LimsMolTableModel molTableModel;
    Vector<Batch> batches = new Vector<Batch>();
    Batch noBatch = new Batch(-1,"Anonymous","Compound Not Assigned.",new Date(System.currentTimeMillis()));
    JXTable batchTable;
    JScrollPane molTableScrollPane;
    Batch currentBatch = noBatch;
    DesignProgressMonitor progressMonitor = new DesignProgressMonitor(LIMSPanel.this,"Progress","Progress",0,100);
    public LIMSPanel() {
        super(new BorderLayout());
        add(buildToolBar(),BorderLayout.NORTH);
        batches.add(noBatch);
        batches.addAll(LimsDAO.getInstance().getAllBatches());
        limsMolecules = new Vector<LimsMolecule>();
        molTableModel = new LimsMolTableModel(limsMolecules);
        molTableModel.addTableModelListener(new TableModelListener() {
            @Override
            public void tableChanged(TableModelEvent e) {
                if(molTableScrollPane!=null){
                    molTableScrollPane.setBorder(new TitledBorder(String.format("%d molecules",limsMolecules.size())));
                }
            }
        });

        molTable2D = new MolTable2D(molTableModel);
        int viewIndex = molTable2D.convertColumnIndexToView(LimsMolecule.PROJECT_IDX);
        molTable2D.getColumn(viewIndex).setCellEditor(new ProjectComboboxEditor(LimsDAO.getInstance().getMyProjectNames()));

        batchTableModel = new BatchTableModel(batches);
        batchTableModel.addTableModelListener(new TableModelListener() {
            @Override
            public void tableChanged(TableModelEvent e) {
                batchTableModel.clearCache();
            }
        });
        JideSplitPane p = new JideSplitPane(JideSplitPane.HORIZONTAL_SPLIT);
        batchTable = new JXTable(batchTableModel);
        batchTable.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        batchTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                int selectedRow = batchTable.getSelectedRow();
                if(selectedRow>=0){
                    int modelRow = batchTable.convertRowIndexToModel(selectedRow);
                    currentBatch = batches.get(modelRow);
                    updateCurrentBatch(currentBatch);
                }
            }
        });
        p.add(new JScrollPane(batchTable));
        molTableScrollPane = new JScrollPane(molTable2D);
        p.add(molTableScrollPane);
        p.setProportionalLayout(true);
        p.setProportions(new double[]{0.3});
        add(p,BorderLayout.CENTER);

        setPreferredSize(new Dimension(1280,1024));
        updateCurrentBatch(currentBatch);

    }

    synchronized void updateCurrentBatch(final Batch batch) {
        SwingWorker sw = new SwingWorker() {
            @Override
            protected Object doInBackground() throws Exception {
                if(batch.getBatch_id()==-1){
                    try {
                        Vector<LimsMolecule> molsNoBatch = LimsDAO.getInstance().getMolsNoBatch();
                        limsMolecules.clear();
                        limsMolecules.addAll(molsNoBatch);
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                    }
                }else{
                    try {
                        Vector<LimsMolecule> molsWithBatch = LimsDAO.getInstance().getMolsByBatch(batch.getBatch_id());
                        limsMolecules.clear();
                        if (molsWithBatch != null) {
                            limsMolecules.addAll(molsWithBatch);
                        }
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                    }
                }
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                    molTableModel.fireTableDataChanged();
                    molTable2D.fixTableFormat();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    e.printStackTrace();
                } finally {
                }
            }
        };
        sw.execute();
    }

    JToolBar buildToolBar(){
        JToolBar toolBar = new JToolBar();
        JButton addNewMolBtn = new JButton("Add New Molecule");
        addNewMolBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                LimsCompoundDialog dialog = new LimsCompoundDialog();
                dialog.setLocationRelativeTo(LIMSPanel.this);
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    if(currentBatch.getBatch_id()==-1) {
                        updateCurrentBatch(currentBatch);
                    }
                    batchTableModel.fireTableDataChanged();
                    JOptionPane.showMessageDialog(LIMSPanel.this,"Molecule added.");
                }
            }
        });
        toolBar.add(addNewMolBtn);

        JButton addBatchBtn = new JButton("Add New Batch");
        addBatchBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                AddBatchDialog dialog = new AddBatchDialog();
                dialog.setLocationRelativeTo(LIMSPanel.this);
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    Batch batch = dialog.getBatch();
                    if(batch!=null){
                        batches.add(batch);
                        batchTableModel.fireTableRowsInserted(batches.size()-1,batches.size()-1);
                    }
                }
            }
        });
        toolBar.add(addBatchBtn);

        JButton deleteBatchBtn = new JButton("Delete Batch");
        deleteBatchBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(currentBatch.getBatch_id()!=-1){
                    try {
                        LimsDAO.getInstance().deleteBatch(currentBatch.getBatch_id());
                        batches.clear();
                        batches.add(noBatch);
                        batches.addAll(LimsDAO.getInstance().getAllBatches());
                        batchTableModel.fireTableDataChanged();
                        limsMolecules.clear();
                        molTableModel.fireTableDataChanged();
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                    }
                }
            }
        });
        toolBar.add(deleteBatchBtn);

        JButton addMolToBatchBtn = new JButton("Add to batch");
        addMolToBatchBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Vector<Integer> selectedMolIds = new Vector<>();
                for(LimsMolecule p:limsMolecules){
                    if(p.getSelected()&&p.getCmpd_id()!=-1){
                        selectedMolIds.add(p.getCmpd_id());
                    }
                }
                if(selectedMolIds.size()==0){
                    JOptionPane.showMessageDialog(LIMSPanel.this,"No molecules selected.");
                    return;
                }
                Vector<Batch> allBatches = LimsDAO.getInstance().getAllBatches();
                if(allBatches.size()==0){
                    JOptionPane.showMessageDialog(LIMSPanel.this,"No batches available.");
                    return;
                }
                Batch[] batchesArray = allBatches.toArray(new Batch[allBatches.size()]);
                Batch batch = (Batch)JOptionPane.showInputDialog(LIMSPanel.this, "Add compounds to batch", "Add to batch", JOptionPane.INFORMATION_MESSAGE, null, batchesArray, batchesArray[allBatches.size()-1]);
                if(batch!=null&&batch.getBatch_id()>=0){
                    try {
                        int row = batchTable.getSelectedRow();
                        if(row<0||batchTable.convertRowIndexToModel(row)<0||batches.get(batchTable.convertRowIndexToModel(row))==null){
                            return;
                        }
                        row = batchTable.convertRowIndexToModel(row);
                        LimsDAO.getInstance().deleteFromBatches(batches.get(row).getBatch_id(),selectedMolIds);
                        LimsDAO.getInstance().addCompoundsToBatch(batch.getBatch_id(),selectedMolIds);
                        updateCurrentBatch(currentBatch);
                        batchTableModel.fireTableDataChanged();
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(LIMSPanel.this,e1.getMessage());
                    }
                }
            }
        });
        toolBar.add(addMolToBatchBtn);


        JButton deleteFromBatchBtn = new JButton("Delete from batch");
        deleteFromBatchBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(currentBatch!=null){
                    Vector<Integer> selectedMoleculeIds = new Vector<>();
                    for(LimsMolecule mol:limsMolecules){
                        if(mol.getSelected()&&mol.getCmpd_id()!=-1){
                            selectedMoleculeIds.add(mol.getCmpd_id());
                        }
                    }
                    if(selectedMoleculeIds.size()==0){
                        return;
                    }
                    if(currentBatch.getBatch_id()==-1){
                        int selection = JOptionPane.showConfirmDialog(LIMSPanel.this, "Molecules will be deleted from database, are you sure?");
                        if(selection!=JOptionPane.YES_OPTION){
                            return;
                        }
                        try {
                            LimsDAO.getInstance().deleteCompounds(selectedMoleculeIds);
                            updateCurrentBatch(currentBatch);
                            batchTableModel.fireTableDataChanged();
                        } catch (SQLException e1) {
                            e1.printStackTrace();
                        }
                        return;
                    }
                    try {
                        LimsDAO.getInstance().deleteFromBatches(currentBatch.getBatch_id(),selectedMoleculeIds);
                        updateCurrentBatch(currentBatch);
                        batchTableModel.fireTableDataChanged();
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                    }
                }
            }
        });
        toolBar.add(deleteFromBatchBtn);

        return toolBar;
    }

    public JMenuBar getMenuBar(){
        JMenuBar menuBar = new JMenuBar();
        JMenu fileMenu = new JMenu("File");
        JMenuItem importItem = new JMenuItem("Import sdf ...");
        importItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser fc = new JFileChooser();
                fc.setFileFilter(new FileFilter() {
                    @Override
                    public boolean accept(File f) {
                        if(f!=null){
                            if(f.isDirectory()){
                                return true;
                            }
                            if(f.getName().endsWith(".sdf")){
                                return true;
                            }
                        }
                        return false;
                    }

                    @Override
                    public String getDescription() {
                        return "SDF file from Electronic notebook.";
                    }
                });
                int option = fc.showOpenDialog(LIMSPanel.this);
                if(option==JFileChooser.APPROVE_OPTION&&fc.getSelectedFile()!=null&&fc.getSelectedFile().canRead()){
                    final DesignProgressMonitor monitor = new DesignProgressMonitor(LIMSPanel.this,"Progress","Progress",0,100, true);
                    monitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    final String sdfFile = fc.getSelectedFile().getAbsolutePath();
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            Vector<LimsMolecule> pmols = new Vector<LimsMolecule>();
                            oemolistream ifs = new oemolistream();
                            ifs.open(sdfFile);
                            OEGraphMol oemol = new OEGraphMol();
                            String scientist = System.getProperty("user.name");
                            Long timestamp = System.currentTimeMillis();
                            while(oechem.OEReadMolecule(ifs,oemol)){
                                if(monitor.isCanceled()){
                                    break;
                                }
                                PropertyMolecule p = new PropertyMolecule(oemol);
                                String compound_id = "";
                                if(oechem.OEHasSDData(oemol,"Compound ID")) {
                                    compound_id = oechem.OEGetSDData(oemol, "Compound ID");
                                }
                                p.addProperty("Compound ID", compound_id);

                                String theory_mass = "";
                                if(oechem.OEHasSDData(oemol,"Theory Mass")) {
                                    theory_mass = oechem.OEGetSDData(oemol, "Theory Mass");
                                }
                                p.addProperty("Theory Mass", theory_mass);

                                String project = "";
                                if(oechem.OEHasSDData(oemol,"Project")){
                                    project = oechem.OEGetSDData(oemol,"Project");
                                }
                                p.addProperty("Exact MW",""+oechem.OECalculateMolecularWeight(p.getMol(),true));

                                pmols.add(new LimsMolecule(p,-1,compound_id,theory_mass,project,scientist,null, null,timestamp));
                            }
                            if(monitor.isCanceled()||pmols.isEmpty()){
                                return null;
                            }else {
                                LimsDAO.getInstance().insertCompounds(pmols, new ProgressReporter() {
                                    @Override
                                    public void reportProgress(String note, int progress) {
                                        Vector v = new Vector();
                                        v.add(note);
                                        v.add(progress);
                                        publish(v);
                                    }
                                });
                                return pmols;
                            }
                        }

                        @Override
                        protected void process(List chunks) {
                            Vector v = (Vector) chunks.get(chunks.size()-1);
                            String note = (String)v.get(0);
                            int progress = (Integer) v.get(1);
                            monitor.setNote(note);
                            monitor.setProgress(progress);

                        }

                        @Override
                        protected void done() {
                            try {
                                Vector<LimsMolecule> pmols = (Vector<LimsMolecule>)get();
                                if(pmols==null){
                                    JOptionPane.showMessageDialog(LIMSPanel.this,"Canceled!");
                                    return;
                                }
                                limsMolecules.clear();
                                limsMolecules.addAll(pmols);
                                batchTableModel.fireTableDataChanged();
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                            } finally {
                                monitor.close();
                            }
                        }
                    };
                    sw.execute();
                }



            }
        });
        fileMenu.add(importItem);
        menuBar.add(fileMenu);

        JMenuItem exportItem = new JMenuItem("Export current batch ...");
        exportItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("SDF","sdf"), new FileFunctor() {
                    @Override
                    public void execute(final File file) {
                        SwingWorker sw = new SwingWorker() {
                            @Override
                            protected Object doInBackground() throws Exception {
                                oemolostream ofs = new oemolostream();
                                ofs.SetFormat(OEFormat.SDF);
                                ofs.open(file.getAbsolutePath());
                                int progress = 0;
                                for(LimsMolecule c:limsMolecules){
                                    Vector v = new Vector();
                                    v.add(String.format("Saving molecule No. %d",progress));
                                    v.add(100*progress/limsMolecules.size());
                                    publish(v);
                                    OEGraphMol mol1 = new OEGraphMol(c.getPropertyMolecule().getMol());
                                    oechem.OEWriteMolecule(ofs,mol1);
                                    progress ++;
                                }
                                ofs.close();
                                return null;
                            }

                            @Override
                            protected void process(List chunks) {
                                Vector v = (Vector) chunks.get(chunks.size()-1);
                                String note = (String)v.get(0);
                                int progress = (Integer) v.get(1);
                                progressMonitor.setNote(note);
                                progressMonitor.setProgress(progress);
                            }

                            @Override
                            protected void done() {
                                try {
                                    get();
                                    progressMonitor.close();
                                    JOptionPane.showMessageDialog(LIMSPanel.this, "SDF file saved.");
                                } catch (Exception e1) {
                                    e1.printStackTrace();
                                    progressMonitor.close();
                                    JOptionPane.showMessageDialog(LIMSPanel.this,e1.getMessage());
                                }
                            }
                        };
                        sw.execute();


                    }
                });
            }
        });
        fileMenu.add(exportItem);

        JMenuItem chiralEnumerationItem = new JMenuItem("Enumerate Chiral centers ...");
        chiralEnumerationItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("SDF","sdf"), new FileFunctor() {
                    @Override
                    public void execute(final File file) {
                        SwingWorker sw = new SwingWorker() {
                            @Override
                            protected Object doInBackground() throws Exception {
                                oemolostream ofs = new oemolostream();
                                ofs.SetFormat(OEFormat.SDF);
                                ofs.open(file.getAbsolutePath());
                                int progress = 0;
                                for(LimsMolecule c:limsMolecules){
                                    Vector v = new Vector();
                                    v.add(String.format("Saving molecule No. %d",progress));
                                    v.add(100*progress/limsMolecules.size());
                                    publish(v);
                                    OEGraphMol mol = c.getPropertyMolecule().getMol();
                                    Vector<OEGraphMol> chiralMols = OEChemFunc.getChiralMols(mol, true);
                                    int n =0;
                                    for(OEGraphMol chiralMol:chiralMols){
                                        oechem.OECopySDData(mol,chiralMol);
                                        chiralMol.SetTitle(String.format("%s_E%d",mol.GetTitle(),++n));
                                        System.out.println(oechem.OEMolToSmiles(chiralMol));
                                    }
                                    for(int i=0;i<4;i++){
                                        int compound_id = i+1;
                                        if(i<chiralMols.size()){
                                            OEGraphMol mol1 = chiralMols.get(i);
                                            mol1 = OEChemFunc.getInstance().getMol2D(mol1);
                                            oechem.OESetSDData(mol1,"Compound_No",""+compound_id);
                                            oechem.OEWriteMolecule(ofs, mol1);
                                        }else{
                                            OEGraphMol tmpMol = new OEGraphMol();
                                            oechem.OESetSDData(tmpMol,"Compound_No",""+compound_id);
                                            oechem.OEWriteMolecule(ofs,tmpMol);
                                            System.out.println(ChemFunc.getMolString(tmpMol));
                                        }
                                    }
                                    progress ++;
                                }
                                ofs.close();
                                return null;
                            }

                            @Override
                            protected void process(List chunks) {
                                Vector v = (Vector) chunks.get(chunks.size()-1);
                                String note = (String)v.get(0);
                                int progress = (Integer) v.get(1);
                                progressMonitor.setNote(note);
                                progressMonitor.setProgress(progress);
                            }

                            @Override
                            protected void done() {
                                try {
                                    get();
                                    progressMonitor.close();
                                    JOptionPane.showMessageDialog(LIMSPanel.this, "SDF file saved.");
                                } catch (Exception e1) {
                                    e1.printStackTrace();
                                    progressMonitor.close();
                                    JOptionPane.showMessageDialog(LIMSPanel.this,e1.getMessage());
                                }
                            }
                        };
                        sw.execute();


                    }
                });            }
        });
        fileMenu.add(chiralEnumerationItem);

        JMenu editMenu = new JMenu("Edit");
        menuBar.add(editMenu);

        JMenuItem selectAllItem = new JMenuItem("Select All");
        selectAllItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(LimsMolecule m:limsMolecules){
                    m.setSelected(true);
                }
                molTableModel.fireTableDataChanged();
            }
        });

        editMenu.add(selectAllItem);

        JMenuItem selectNoneItem = new JMenuItem("Deselect All");
        selectNoneItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(LimsMolecule m:limsMolecules){
                    m.setSelected(false);
                }
                molTableModel.fireTableDataChanged();
            }
        });

        editMenu.add(selectNoneItem);

        JMenuItem inverseSelectionItem = new JMenuItem("Invert Selection");
        inverseSelectionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(LimsMolecule m:limsMolecules){
                    if(m.getSelected()){
                        m.setSelected(false);
                    }else{
                        m.setSelected(true);
                    }
                }
                molTableModel.fireTableDataChanged();
            }
        });

        editMenu.add(inverseSelectionItem);

        JMenuItem selectHighlightItem = new JMenuItem("Select hightlighted");
        selectHighlightItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(Integer row:molTable2D.getSelectedRows()){
                    if(row>=0){
                        row = molTable2D.convertRowIndexToModel(row);
                        limsMolecules.get(row).setSelected(true);
                    }
                }
                molTableModel.fireTableDataChanged();
            }
        });

        editMenu.add(selectHighlightItem);

        return menuBar;
    }


    class ProjectComboboxEditor extends AbstractCellEditor implements TableCellEditor,ActionListener{
        String project;
        Vector<String> all_projects;
        int row = -1;

        public ProjectComboboxEditor(Vector<String> all_projects) {
            this.all_projects = all_projects;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            JComboBox<String> projectCB = (JComboBox<String>) e.getSource();
            this.project = projectCB.getItemAt(projectCB.getSelectedIndex());
            fireEditingStopped();
            molTableModel.fireTableRowsUpdated(row,row);
        }

        @Override
        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
            if (value instanceof String) {
                this.project = (String)value;
            }
            this.row = row;
            JComboBox<String> projectCB = new JComboBox<>(all_projects);
            projectCB.setSelectedItem(value);
            projectCB.addActionListener(this);

//            if (isSelected) {
//                statusCB.setBackground(table.getSelectionBackground());
//            } else {
//                statusCB.setBackground(table.getSelectionForeground());
//            }

            return projectCB;
        }

        @Override
        public boolean isCellEditable(EventObject e) {
            if (e instanceof MouseEvent) {
                return ((MouseEvent)e).getClickCount() >= 2;
            }
            return true;
        }

        @Override
        public Object getCellEditorValue() {
            return this.project==null?"":this.project;
        }
    }




    public static void main(String[] args) {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        LIMSPanel p = new LIMSPanel();
        frame.getContentPane().add(p);
        frame.setJMenuBar(p.getMenuBar());
        frame.pack();
        frame.setVisible(true);

    }


}
