package com.insilico.application.insilicotools.gui.widget;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.license.LicenseManager;
import chemaxon.license.LicenseProcessingException;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.Core;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.gui.table.MyCoreTableModel;
import com.insilico.application.insilicotools.gui.table.MolTable2D;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.gui.util.TableUtils;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemWebLicenseInstaller;
import com.insilico.application.insilicotools.gui.*;
import com.jidesoft.swing.JideScrollPane;
import com.jidesoft.swing.JideSplitPane;
import com.jidesoft.swing.PartialEtchedBorder;
import openeye.oechem.*;
import org.jdesktop.swingx.JXImageView;
import org.jdesktop.swingx.JXTable;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileFilter;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 2/13/17.
 */
public class MyCoreBrowsingPanel extends JPanel {
    MyCoreTableModel tableModel;
    JXTable table;
    Vector<Core> myCores;
    JXImageView viewPanel;
    MolTable2D table2D;
    Vector<PropertyMolecule> exampleMols;
    PropertyMolTableModel molTableModel;
    final static int max_num_of_examples = 100;
    JButton addExamplesBtn ;
    JButton removeBtn;
    JButton editBtn;
    JButton exportSdfBtn;
    Vector<String> properties;
    CompoundInputDialog compoundInputDialog;
    String smarts;
    DesignProgressMonitor progressMonitor = new DesignProgressMonitor(MyCoreBrowsingPanel.this,"Progress","Progress",0,100);


    public JMenuBar getJMenuBar(){
        JMenuBar menuBar = new JMenuBar();
        JMenu taskMenu = new JMenu("Task");
        JMenuItem batchItem = new JMenuItem("Batch loading cores ...");
        batchItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MyCoreBatchRegistrationDialog dialog = new MyCoreBatchRegistrationDialog();
                dialog.setLocationRelativeTo(MyCoreBrowsingPanel.this);
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    final String libraryNameTag = dialog.getLibraryNameTag();
                    final String sdfPathName = dialog.getSdfFilePath();
                    final String source = dialog.getSource();
                    progressMonitor.setMillisToDecideToPopup(0);
                    progressMonitor.setMillisToPopup(0);
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            FrontierDAO.getInstance().insertMyCoresBatchWrapper(sdfPathName,"chemist",source,libraryNameTag);
                            return null;
                        }

                        @Override
                        protected void done() {
                            try {
                                get();
                                reload();
                                progressMonitor.close();
                                JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,"Libraries imported.");
                            } catch (Exception e1) {
                                progressMonitor.close();
                                e1.printStackTrace();
                            }
                        }
                    };
                    sw.execute();
                }
            }
        });
        taskMenu.add(batchItem);
        menuBar.add(taskMenu);
        return menuBar;
    }

    public MyCoreBrowsingPanel(final Vector<Core> myCores) {
        super(new BorderLayout());
        JideSplitPane p1 = new JideSplitPane(JideSplitPane.HORIZONTAL_SPLIT);
        compoundInputDialog = new CompoundInputDialog();
        this.myCores = myCores;
        tableModel = new MyCoreTableModel(myCores);
        viewPanel = new JXImageView();

        exampleMols = new Vector<PropertyMolecule>();
        properties = new Vector<String>();
        properties.add("molecular weight");
        properties.add("2d PSA");
        properties.add("CLogP");

        molTableModel = new PropertyMolTableModel(exampleMols, properties);
        table2D = new MolTable2D(molTableModel);

        buildTable();

        JideScrollPane scrollPane = new JideScrollPane(table);
        final JPanel tablePanel = new JPanel(new BorderLayout());
        tablePanel.add(scrollPane,BorderLayout.CENTER);
        tablePanel.setBorder(new PartialEtchedBorder());
        JPanel btnPanel1 = new JPanel();
        final JCheckBox smartsCB = new JCheckBox("Substructure",false);
        smartsCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(smartsCB.isSelected()){
                    compoundInputDialog.setLocationRelativeTo(MyCoreBrowsingPanel.this);
                    compoundInputDialog.setVisible(true);
                    if(compoundInputDialog.isCommitted()&&(!compoundInputDialog.getMolecule().isEmpty())){
                        Molecule mol = compoundInputDialog.getMolecule();
                        try {
                            smarts = MolExporter.exportToFormat(mol, "smarts:u,a");
                        } catch (IOException e1) {
                            smartsCB.setSelected(false);
                            smarts = null;
                            e1.printStackTrace();
                        }
                        try {
                            reload();
                        } catch (SQLException e1) {
                            e1.printStackTrace();
                        }
                    }else{
                        smartsCB.setSelected(false);
                        smarts = null;
                    }
                    System.out.println("Selected");
                }else{
                    smarts = null;
                    try {
                        reload();
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                    }
                }
            }
        });
        btnPanel1.add(smartsCB);


        tablePanel.add(btnPanel1,BorderLayout.SOUTH);
        p1.setProportionalLayout(true);
        p1.add(tablePanel);
        JPanel p3 = new JPanel(new BorderLayout());
        p3.add(new JideScrollPane(table2D),BorderLayout.CENTER);
        JPanel btnPanel2 = new JPanel();
        addExamplesBtn = new JButton("Add sample molecules");
        addExamplesBtn.setEnabled(false);
        addExamplesBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser fc = new JFileChooser(InSlilicoPanel.getInstance().getCurrentDirectory());
                fc.setFileFilter(new FileFilter() {
                    @Override
                    public boolean accept(File f) {
                        if(f.getName().endsWith(".sdf")){
                            return true;
                        }
                        if(f.isDirectory()){
                            return true;
                        }
                        return false;
                    }

                    @Override
                    public String getDescription() {
                        return "SDF file";
                    }
                });
                if(fc.showOpenDialog(MyCoreBrowsingPanel.this)==JFileChooser.APPROVE_OPTION){
                    oemolistream ifs = new oemolistream();
                    File selectedFile = fc.getSelectedFile();
                    if(selectedFile!=null&&selectedFile.canRead()) {
                        ifs.open(selectedFile.getAbsolutePath());
                        OEGraphMol oemol = new OEGraphMol();
                        exampleMols.clear();
                        int n = 0;
                        while (oechem.OEReadMolecule(ifs, oemol)) {
                            PropertyMolecule pmol = new PropertyMolecule(oemol);
                            exampleMols.add(pmol);
                            n += 1;
                            if(n>max_num_of_examples){
                                break;
                            }
                        }
                        ChemFunc.calculateOEProperty(exampleMols);
                        for(Core core:myCores) {
                            if(core.getSelected()) {
                                try {
                                    FrontierDAO.getInstance().insertSampleMolecules(core.getId(),exampleMols,null);
                                } catch (SQLException e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this, e1.getMessage());
                                }
                            }
                        }
                        molTableModel.setSelectedProperties(properties);
                        molTableModel.fireTableStructureChanged();
                    }
                }
            }
        });
        btnPanel2.add(addExamplesBtn);
        p3.add(btnPanel2,BorderLayout.SOUTH);

        JideSplitPane splitPane = new JideSplitPane(JideSplitPane.VERTICAL_SPLIT);
        splitPane.add(p3);
        splitPane.add(viewPanel);
        splitPane.setProportionalLayout(true);
        splitPane.setProportions(new double[]{0.5});

        p1.add(splitPane);
        p1.setProportions(new double[]{0.5});


        JPanel btnPanel = new JPanel();
        JButton addBtn = new JButton("Add New Library");
        addBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                final JDialog dialog = new JDialog();
                final MyCoreInputPanel inputPanel = new MyCoreInputPanel();
                inputPanel.addPropertyChangeListener(new PropertyChangeListener() {
                    @Override
                    public void propertyChange(PropertyChangeEvent evt) {
                        if(evt.getPropertyName().equals("LibraryUpdated")){
                            dialog.setVisible(false);
                            Core core = inputPanel.getCore();
                            myCores.insertElementAt(core,0);
                            tableModel.fireTableDataChanged();
                        }
                    }
                });
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.setLocationRelativeTo(MyCoreBrowsingPanel.this);
                dialog.setModal(true);
                dialog.setSize(new Dimension(600,700));
                dialog.getContentPane().add(inputPanel);
                dialog.setVisible(true);

            }
        });
        btnPanel.add(addBtn);

        editBtn = new JButton("Edit selected");
        editBtn.setEnabled(false);
        editBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(myCores==null||myCores.isEmpty()){
                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,"No core available.");
                    return;
                }
                Core selectedCore = null;
                for(Core core:myCores){
                    if(core.getSelected()){
                        selectedCore = core;
                        break;
                    }
                }
                if(selectedCore==null){
                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,"No core selected.");
                    return;
                }
                final JDialog dialog = new JDialog();
                MyCoreInputPanel inputPanel = new MyCoreInputPanel(selectedCore);
                inputPanel.addPropertyChangeListener(new PropertyChangeListener() {
                    @Override
                    public void propertyChange(PropertyChangeEvent evt) {
                        if(evt.getPropertyName().equals("LibraryUpdated")){
                            dialog.setVisible(false);
                            tableModel.fireTableDataChanged();
                        }
                    }
                });
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.setLocationRelativeTo(MyCoreBrowsingPanel.this);
                dialog.setModal(true);
                dialog.setSize(new Dimension(600,700));
                dialog.getContentPane().add(inputPanel);
                dialog.setVisible(true);
            }
        });
        btnPanel.add(editBtn);

        removeBtn = new JButton("Remove selected");
        removeBtn.setEnabled(false);
        removeBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(myCores.isEmpty()){
                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,"No library available.");
                }
                int decision = JOptionPane.showConfirmDialog(MyCoreBrowsingPanel.this,"Are you sure?");
                if(decision!=JOptionPane.YES_OPTION){
                    return;
                }
                int n = 0;
                Vector<Core> coresToRemove = new Vector<Core>();
                for(Core core:myCores){
                    if(core.getSelected()){
                        try {
                            FrontierDAO.getInstance().removeMyCore(core.getId());
                            coresToRemove.add(core);
                            n += 1;
                        } catch (SQLException e1) {
                            e1.printStackTrace();
                        }
                    }
                }
                if(!coresToRemove.isEmpty()){
                    for(Core core:coresToRemove){
                        myCores.remove(core);
                    }
                }
                tableModel.fireTableDataChanged();
                if(n>1) {
                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this, String.format("%d libraries removed.", n));
                }else{
                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this, String.format("%d library removed.", n));
                }
                if(n==1){
                    exampleMols.clear();
                    molTableModel.fireTableStructureChanged();
                }

            }
        });
        btnPanel.add(removeBtn);

        exportSdfBtn = new JButton("Export SDF ...");
        exportSdfBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                saveIdeaAsSdf();
            }
        });
        btnPanel.add(exportSdfBtn);

        add(p1,BorderLayout.CENTER);
        add(btnPanel,BorderLayout.SOUTH);
    }

    void saveIdeaAsSdf() {
        if(myCores.size()==0){
            JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,"No cores are available.");
            return;
        }
        FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("SDF file","sdf"), new FileFunctor() {
            @Override
            public void execute(final File file) {
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        oemolostream ofs = new oemolostream();
                        ofs.SetFormat(OEFormat.SDF);
                        ofs.open(file.getAbsolutePath());
                        int progress = 0;
                        for(Core c: myCores){
                            Vector v = new Vector();
                            v.add(String.format("Saving Core No. %d",progress));
                            v.add(100*progress/ myCores.size());
                            publish(v);
                            OEGraphMol mol1 = new OEGraphMol(c.getCore_mol().getMol());
                            for(int idx=2;idx<tableModel.getColumnCount();idx++){
                                String colName = tableModel.getColumnName(idx);
                                Object obj = tableModel.getValueAt(progress, idx);
                                if(obj!=null) {
                                    String value = obj.toString();
                                    if (value != null) {
                                        oechem.OESetSDData(mol1, colName, value);
                                    } else {
                                        oechem.OESetSDData(mol1, colName, "");
                                    }
                                }
                            }
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
                            JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this, "SDF file saved.");
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,e1.getMessage());
                        }finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();

            }
        });
    }

    private void buildTable() {
        table = new JXTable(tableModel);
        table.setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer());
        table.setShowGrid(true);
        table.setGridColor(Color.GRAY);
        table.setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer());
        table.setSortOrderCycle(SortOrder.ASCENDING, SortOrder.DESCENDING, SortOrder.UNSORTED);
        table.setColumnControlVisible(true);
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        table.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        final int selectColIdx = table.convertColumnIndexToModel(0);
        int molColIdx = table.convertColumnIndexToModel(1);
        table.getColumnModel().getColumn(selectColIdx).setMaxWidth(20);

        table.getColumnModel().addColumnModelListener(new TableColumnModelListener() {
            @Override
            public void columnAdded(TableColumnModelEvent e) {

            }

            @Override
            public void columnRemoved(TableColumnModelEvent e) {

            }

            @Override
            public void columnMoved(TableColumnModelEvent e) {

            }

            @Override
            public void columnMarginChanged(ChangeEvent e) {
                if (table.getColumnModel().getColumnCount() > 0) {
                    int molColIdx = table.convertColumnIndexToModel(1);
                    int molTableRowHeight = table.getColumnModel().getColumn(molColIdx).getWidth();
                    table.setRowHeight(molTableRowHeight);
                }
            }

            @Override
            public void columnSelectionChanged(ListSelectionEvent e) {

            }
        });
        table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(!e.getValueIsAdjusting()) {
                    for(Core core: myCores){
                        core.setSelected(false);
                    }
                    int selectedRow = table.getSelectedRow();
                    if(selectedRow>=0&&selectedRow< myCores.size()) {
                        int index = table.convertRowIndexToModel(selectedRow);
                        if (index >= 0) {
                            Core core = tableModel.getCore(index);
                            if (core != null) {
                                core.setSelected(true);
//                                Vector<String> super_users = FrontierDAO.getInstance().getSuper_users();
                                String username = System.getProperty("user.name");
                                if(core.getChemist().equals(username)){
                                    removeBtn.setEnabled(true);
                                    editBtn.setEnabled(true);
                                    addExamplesBtn.setEnabled(true);
                                }else{
                                    removeBtn.setEnabled(false);
                                    editBtn.setEnabled(false);
                                    addExamplesBtn.setEnabled(false);
                                }
                                byte[] cdx = core.getCdx();
                                try {
                                    Molecule mol = MolImporter.importMol(cdx, "cdx", "utf-8");
                                    byte[] png = MolExporter.exportToBinFormat(mol, "png:w600,h400");
                                    viewPanel.setImage(ImageIO.read(new ByteArrayInputStream(png)));
                                } catch (IOException e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this, e1.getMessage());
                                }
                                try {
                                    exampleMols.clear();
                                    Vector<PropertyMolecule> mols = FrontierDAO.getInstance().getSampleMolecules(core.getId());
                                    if(!mols.isEmpty()) {
                                        ChemFunc.calculateOEProperty(mols);
                                    }
                                    exampleMols.addAll(mols);
                                    molTableModel.setSelectedProperties(properties);
                                    molTableModel.fireTableStructureChanged();
                                } catch (SQLException e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(MyCoreBrowsingPanel.this,e1.getMessage());
                                }
                            }
                        }
                    }
                    table.repaint();
                }
            }
        });
    }

    void reload() throws SQLException{
        final int selectedRow = table.getSelectedRow();
        myCores.clear();
        SwingWorker sw = new SwingWorker() {
            @Override
            protected Object doInBackground() throws Exception {
                if(smarts==null) {
                    return FrontierDAO.getInstance().getAllMyCores();
                }else{
                    return FrontierDAO.getInstance().getMyCoresBySubstructure(smarts);
                }
            }

            @Override
            protected void done() {
                try{
                    Vector<Core> cores = (Vector<Core>)get();
                    myCores.clear();
                    myCores.addAll(cores);
                    tableModel.fireTableDataChanged();
                    if(selectedRow>=0) {
                        TableUtils.scrollToVisible(table, selectedRow);
                    }

                } catch (InterruptedException e) {
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    e.printStackTrace();
                }
            }
        };
        sw.execute();
    }

    public static void main(String[] args) {
        try {
            OEChemWebLicenseInstaller.loadOELicenseFromWeb();
//                    LicenseManager.setLicenseFile("/Users/jfeng1/.chemaxon/license.cxl");
            LicenseManager.setLicenseFile("http://10.74.2.128:8080/insilicotools/license.cxl");
        } catch (IOException e) {
            JOptionPane.showMessageDialog(null,e.getMessage());
            return;
        } catch (LicenseProcessingException e) {
            JOptionPane.showMessageDialog(null,e.getMessage());
            return;
        }

        JFrame frame = new JFrame("test");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(new Dimension(1024,768));
        Vector<Core> allMyCores = null;
        try {
            allMyCores = FrontierDAO.getInstance().getAllMyCores();
        } catch (SQLException e) {
            e.printStackTrace();
            return;
        }
        MyCoreBrowsingPanel p = new MyCoreBrowsingPanel(allMyCores);
        frame.setJMenuBar(p.getJMenuBar());
        frame.getContentPane().add(p);
        frame.setVisible(true);
    }
}
