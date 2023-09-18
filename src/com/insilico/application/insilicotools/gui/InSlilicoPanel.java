package com.insilico.application.insilicotools.gui;

import chemaxon.formats.MolExporter;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.InSilicoToolOptions;
import com.insilico.application.insilicotools.chart.ClearancePlot;
import com.insilico.application.insilicotools.chart.EganEgg;
import com.insilico.application.insilicotools.chart.PrincipalMomentOfInertiaPlot;
import com.insilico.application.insilicotools.cmdline.BatchDescriptorGenerationPanel;
import com.insilico.application.insilicotools.data.*;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.insilico.application.insilicotools.gui.dialog.MultiChoicesDialog;
import com.insilico.application.insilicotools.gui.dialog.MultiInputDialog;
import com.insilico.application.insilicotools.gui.dialog.MultiInputTextDialog;
import com.insilico.application.insilicotools.gui.filter.FilterDialog;
import com.insilico.application.insilicotools.gui.lims.LIMSPanel;
import com.insilico.application.insilicotools.gui.modeling.ConfSearchInputPanel;
import com.insilico.application.insilicotools.gui.modeling.MolDockingWizard;
import com.insilico.application.insilicotools.gui.modeling.MolOverlayWizard;
import com.insilico.application.insilicotools.gui.table.MatrixCellTable;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.gui.util.ImageUtil;
import com.insilico.application.insilicotools.gui.widget.*;
import com.insilico.application.insilicotools.inSilicoTools;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.google.common.base.Joiner;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import com.jidesoft.range.IntegerRange;
import com.jidesoft.range.NumericRange;
import com.jidesoft.range.Range;
import com.jidesoft.swing.JideScrollPane;
import com.jidesoft.swing.JideTabbedPane;
import openeye.oechem.*;
//import openeye.oedepict.oedepict;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.jdesktop.swingx.JXStatusBar;
import org.jdesktop.swingx.JXTable;
import org.jdesktop.swingx.search.TableSearchable;
import com.insilico.application.insilicotools.util.ChemFunc;
import org.jfree.ui.RefineryUtilities;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.*;
import java.net.*;
import java.sql.SQLException;
import java.util.*;
import java.util.List;
import java.util.Timer;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 9/10/15.
 */
public class InSlilicoPanel extends JPanel {
    PropertyMolTableModel tableModel;
    Vector<PropertyMolecule> molecules;
    File currentDirectory;
    DesignProgressMonitor progressMonitor;
    JXTable molTable;
    EganEgg plot;
    PrincipalMomentOfInertiaPlot pmiPlot;
    ClearancePlot clearancePlot;
    PropertySelectionDialog propertySelectionDialog;
    Vector<String> existingTags;
    JXStatusBar statusBar;
    JLabel numMolLabel;
    MatrixMolTableModel matrixMolTableModel;
    MatrixCellTable matrixTable;
    JComboBox columnCountCB;
    JSlider columnWidthSlider;
    JButton saveAsPdfBtn;
    CompoundInputDialog sketcher;
    public static InSlilicoPanel _this;
    int molTableRowHeight = 150;
    JideTabbedPane tabbedPane;
    JPanel workspacePanel;
    JideTabbedPane workspaceTabbedPanel;
    LIMSPanel limsPanel;
    MyCoreBrowsingPanel coresPanel;
    public final static String ENAMINE_PREFERRED = "Enamine Preferred";
    public final static String ENAMINE = "Enamine";
    public final static String FRONTIER = "Frontier";
    public final static String PHARMARON = "Pharmaron";
    public final static String MARKET_SELECT = "Market Select";
    public final static String SKYHAWK = "Cellarity";

    public static InSlilicoPanel getInstance(){
        if(_this==null){
            _this = new InSlilicoPanel();
        }
        return _this;
    }

    public String getUserName(){
        return System.getProperty("user.name");
    }

    public File getCurrentDirectory() {
        return currentDirectory;
    }

    public void setCurrentDirectory(File currentDirectory) {
        this.currentDirectory = currentDirectory;
    }

    private InSlilicoPanel() {
        super(new BorderLayout());
        molecules = new Vector<PropertyMolecule>();
        progressMonitor = new DesignProgressMonitor(this,"Progress","Loading",0,100);
        existingTags = new Vector<String>();
        sketcher = new CompoundInputDialog();
        workspacePanel = new JPanel(new BorderLayout());
        workspaceTabbedPanel = new JideTabbedPane();
        buildTablePanel();
        buildToolBar();
        buildStatusBar();
        workspacePanel.add(workspaceTabbedPanel,BorderLayout.CENTER);
        tabbedPane.add("Workspace 2D",workspacePanel);
        tabbedPane.setTabClosableAt(0,false);
        tabbedPane.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                JTabbedPane sourceTabbedPane = (JTabbedPane) e.getSource();
                int index = sourceTabbedPane.getSelectedIndex();
                String title = sourceTabbedPane.getTitleAt(index);
                if(title.equals("Workspace 2D")){
                    inSilicoTools.getInstance().setJMenuBar(InSlilicoPanel.getInstance().getMenuBar());
                }else if(title.startsWith("Project")){
                    MyCompoundIdealPanel p = (MyCompoundIdealPanel)sourceTabbedPane.getSelectedComponent();
                    inSilicoTools.getInstance().setJMenuBar(p.getJMenuBar());
                }else if(title.startsWith("Purification")){
                    inSilicoTools.getInstance().setJMenuBar(limsPanel.getMenuBar());
                }else if(title.equalsIgnoreCase("Cores")){
                    inSilicoTools.getInstance().setJMenuBar(coresPanel.getJMenuBar());
                }else{
                    inSilicoTools.getInstance().setJMenuBar(InSlilicoPanel.getInstance().getMenuBar());
                }
            }
        });

    }

    public CompoundInputDialog getSketcher() {
        return sketcher;
    }

    private void buildPopupMenu(){

        molTable.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if(SwingUtilities.isRightMouseButton(e)){
                    JTable source = (JTable)e.getSource();
                    final int row = source.rowAtPoint( e.getPoint() );
                    final int column = source.columnAtPoint( e.getPoint() );
                    if(row<0||column<0){
                        return;
                    }
                    final int modelRow = molTable.convertRowIndexToModel(row);
                    final int modelColumn = molTable.convertColumnIndexToModel(column);
                    final JPopupMenu popupMenu = new JPopupMenu();
                    JMenuItem editMolItem = new JMenuItem("Edit Current Molecule ...");
                    editMolItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(modelRow>=0){
                                final PropertyMolecule propertyMolecule = molecules.get(modelRow);
                                if(propertyMolecule!=null){
                                    sketcher.setMolecule(propertyMolecule);
                                    sketcher.setVisible(true);
                                    if(sketcher.isCommitted()){
                                        Molecule molecule = sketcher.getMolecule();
                                        if(!molecule.isEmpty()){
                                            OEGraphMol oeGraphMol = OEChemFunc.getInstance().convertChemAxonMol(molecule);
                                            propertyMolecule.setOEMol(oeGraphMol);
                                            ChemFunc.calculateOEProperty(new Vector<PropertyMolecule>(){
                                                {
                                                    add(propertyMolecule);
                                                }
                                            });
                                            tableModel.fireTableDataChanged();
                                            matrixMolTableModel.fireTableDataChanged();
                                        }
                                    }
                                }
                            }
                        }
                    });
                    final Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                    JMenuItem copyAsMolItem = new JMenuItem("Copy Curren Molecule as SDF ...");
                    copyAsMolItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(modelRow>=0){
                                PropertyMolecule propertyMolecule = molecules.get(modelRow);
                                String molString = ChemFunc.getMolString(propertyMolecule.getMol());
                                if(molString!=null) {
                                    StringSelection selection = new StringSelection(molString);
                                    clipboard.setContents(selection, null);
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"Molecule copied to clipboard.");
                                    return;
                                }
                            }
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "Failed to copy to clipboard, molecule is not available.");
                        }
                    });
                    JMenuItem copyAsMol3DItem = new JMenuItem("Copy Curren Molecule as 3D SDF ...");
                    copyAsMol3DItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(modelRow>=0){
                                PropertyMolecule propertyMolecule = molecules.get(modelRow);
                                String molString = ChemFunc.getMolString(propertyMolecule.getMol3d());
                                if(molString!=null) {
                                    StringSelection selection = new StringSelection(molString);
                                    clipboard.setContents(selection, null);
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"Molecule copied to clipboard.");
                                    return;
                                }
                            }
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "Failed to copy to clipboard, molecule is not available.");
                        }
                    });
                    JMenuItem copyAsSmilesItem = new JMenuItem("Copy Current Molecule as Smiles ...");
                    copyAsSmilesItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(modelRow>=0){
                                PropertyMolecule propertyMolecule = molecules.get(modelRow);
                                String molName = propertyMolecule.getName()==null?"":" "+propertyMolecule.getName();
                                String molString = propertyMolecule.getSmiles()+molName;
                                if(molString!=null) {
                                    StringSelection selection = new StringSelection(molString);
                                    clipboard.setContents(selection, null);
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"Smiles copied to clipboard.");
                                    return;
                                }
                            }
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "Failed to copy to clipboard, molecule is not available.");
                        }
                    });
                    final int numSelected = getSelectedMolecules().size();
                    JMenuItem copyNameItem = new JMenuItem("Copy Current Molecule Names ...");
                    copyNameItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(modelRow>=0){
                                PropertyMolecule propertyMolecule = molecules.get(modelRow);
                                String molString = propertyMolecule.getName();
                                if(molString!=null) {
                                    StringSelection selection = new StringSelection(molString);
                                    clipboard.setContents(selection, null);
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"Molecule Name copied to clipboard.");
                                    return;
                                }
                            }
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "Failed to copy to clipboard, molecule is not available.");
                        }
                    });
                    JMenuItem copySelectedAsSmilesItem = new JMenuItem("Copy Selected Molecules as Smiles ...");
                    copySelectedAsSmilesItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(numSelected==0){
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule selected.");
                                return;
                            }
                            StringBuilder sb = new StringBuilder();
                            for(PropertyMolecule mol:getSelectedMolecules()){
                                String molName = mol.getName()==null?"":" "+mol.getName();
                                sb.append(mol.getSmiles()+molName + "\n");
                            }
                            StringSelection selection = new StringSelection(sb.toString());
                            clipboard.setContents(selection,null);
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"Smiles copied to clipboard.");
                        }
                    });
                    JMenuItem copySelectedMolNamesItem = new JMenuItem("Copy Selected Molecule Names ...");
                    copySelectedMolNamesItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            if(numSelected==0){
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule selected.");
                                return;
                            }
                            StringBuilder sb = new StringBuilder();
                            for(PropertyMolecule mol:getSelectedMolecules()){
                                sb.append(mol.getName()+"\n");
                            }
                            StringSelection selection = new StringSelection(sb.toString());
                            clipboard.setContents(selection,null);
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"Molecule names copied to clipboard.");
                        }
                    });

                    JMenuItem copySelectedUniqMolNamesItem = new JMenuItem("Copy Unique Selected Molecule Names ...");
                    copySelectedUniqMolNamesItem.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            HashMap<String,Boolean> uniqDict = new HashMap<>();
                            if(numSelected==0){
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule selected.");
                                return;
                            }
                            StringBuilder sb = new StringBuilder();
                            for(PropertyMolecule mol:getSelectedMolecules()){
                                String molName = mol.getName();
                                if(!uniqDict.containsKey(molName)){
                                    uniqDict.put(molName,true);
                                    sb.append(molName + "\n");
                                }
                            }
                            StringSelection selection = new StringSelection(sb.toString());
                            clipboard.setContents(selection,null);
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"Molecule names copied to clipboard.");
                        }
                    });


                    if(numSelected==0){
                        copySelectedAsSmilesItem.setEnabled(false);
                        copySelectedMolNamesItem.setEnabled(false);
                    }
                    popupMenu.add(editMolItem);
                    popupMenu.add(copyAsMolItem);
                    popupMenu.add(copyAsMol3DItem);
                    popupMenu.add(copyAsSmilesItem);
                    popupMenu.add(copyNameItem);
                    popupMenu.add(copySelectedAsSmilesItem);
                    popupMenu.add(copySelectedMolNamesItem);
                    popupMenu.add(copySelectedUniqMolNamesItem);
                    if(modelColumn==1){
                        popupMenu.show(source,e.getX(),e.getY());
                    }
                }
            }
        });


    }

    private void buildStatusBar(){
        statusBar = new JXStatusBar();
        statusBar.setPreferredSize(new Dimension(1000,40));
        numMolLabel = new JLabel("Number of Molecules (0)");
        statusBar.add(numMolLabel);

        JButton alignBtn = new JButton("Align Molecule ...");
        alignBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules==null||molecules.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                    return;
                }
                sketcher.setLocationRelativeTo(InSlilicoPanel.this);
                sketcher.setVisible(true);
                if(sketcher.isCommitted()){
                    progressMonitor.setNote("Aligning molecules ...");
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            String qmolString = MolExporter.exportToFormat(sketcher.getMolecule(), "sdf");
                            OEGraphMol qmol = ChemFunc.getMolFromMolString(qmolString,OEFormat.SDF);
                            OESubSearch subSearch = new OESubSearch();
                            subSearch.Init(qmol,OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
                            for(PropertyMolecule mol:molecules){
                                int progress = 100*molecules.indexOf(mol)/molecules.size();
                                publish(new Integer(progress));
                                oechem.OEPrepareSearch(mol.getMol(),subSearch);
//                              oedepict.OEPrepareAlignedDepiction(mol.getMol(),subSearch);
                                mol.update2DDepiction();
                            }
                            return null;
                        }

                        @Override
                        protected void process(List chunks) {
                            Integer progress = (Integer)chunks.get(chunks.size()-1);
                            progressMonitor.setProgress(progress);
                        }

                        @Override
                        protected void done() {
                            try {
                                get();
                                tableModel.fireTableDataChanged();
                                tableModel.fireTableStructureChanged();
                                fixTableFormat();
                                matrixMolTableModel.fireTableDataChanged();
                                matrixMolTableModel.fireTableStructureChanged();
                                fixMatrixTableFormat();
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            } finally{
                                progressMonitor.setNote("Progress");
                                progressMonitor.close();
                            }
                        }
                    };
                    sw.execute();

                }
            }
        });
        statusBar.add(alignBtn);

        final JLabel label = new JLabel();
        statusBar.add(label);
        Timer timer = new Timer();
        timer.scheduleAtFixedRate(new TimerTask() {
            @Override
            public void run() {
                long freeMemory = Runtime.getRuntime().freeMemory();
                long totalMemory = Runtime.getRuntime().totalMemory();
                long maxMemory = Runtime.getRuntime().maxMemory();
                label.setText(String.format("%d/%d M ",(totalMemory-freeMemory)/1024000, maxMemory/1024000));
            }
        },0,10000);

        JButton gcBtn = new JButton("Clean");
        gcBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                System.gc();
                long freeMemory = Runtime.getRuntime().freeMemory();
                long totalMemory = Runtime.getRuntime().totalMemory();
                long maxMemory = Runtime.getRuntime().maxMemory();
                label.setText(String.format("%d/%d M ",(totalMemory-freeMemory)/1024000, maxMemory/1024000));
            }
        });
        statusBar.add(gcBtn);

        workspacePanel.add(statusBar, BorderLayout.SOUTH);
    }

    private void buildToolBar(){
        JToolBar toolBar = new JToolBar();
        JButton sketchMolBtn = new JButton("Sketch Molecule", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("sketch-icon.png"))));
        sketchMolBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                loadMolFromSketcher();
            }
        });
        toolBar.add(sketchMolBtn);

        JButton fileBtn = new JButton("Load SDF",ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("File-icon_small.png"))));
        fileBtn.addActionListener(new LoadSDFListener());
        toolBar.add(fileBtn);

        JButton pasteBtn = new JButton("Paste Smiles/CY_NUMBER", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("File-icon_small.png"))));
        pasteBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                PasteSmilesDialog dialog = new PasteSmilesDialog();
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setVisible(true);
                if(dialog.getDialogResult()==StandardDialog.RESULT_AFFIRMED){
                    Vector<PropertyMolecule> tmpMols = dialog.getMolecules();
                    if(tmpMols.size()>0){
                        if(molecules.size()>0) {
                            int option = JOptionPane.showConfirmDialog(InSlilicoPanel.this, String.format("%d molecules loaded, overwrite existing molecules?", tmpMols.size()));
                            if (option == JOptionPane.YES_OPTION) {
                                molecules.clear();
                                existingTags.clear();
                            }
                        }
                        for(PropertyMolecule mol:tmpMols){
                            molecules.add(mol);
                        }
                        tableModel.fireTableStructureChanged();
                        tableModel.fireTableDataChanged();
                        fixTableFormat();
                    }
                    updateStatusBar();
                }
            }
        });
        toolBar.add(pasteBtn);

        JButton eggBtn = new JButton("Egan Egg", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("software-egg-icon_small.png"))));
        eggBtn.addActionListener(new drawEganEggListener());
        toolBar.add(eggBtn);

        JButton plotBtn = new JButton("Scatter Plot", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("chart-scatter-icon_small.png"))));
//        toolBar.add(plotBtn);

        JButton excelBtn = new JButton("Export to Excel", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("Excel-icon_small.png"))));
        excelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                    return;
                }
                int option = JOptionPane.showConfirmDialog(InSlilicoPanel.this, "Save mol image?");
                final boolean saveAsImage = option == JOptionPane.YES_OPTION;
                FileUtil.saveToFile(currentDirectory, new FileNameExtensionFilter("Excel file","xls"), new FileFunctor() {
                    @Override
                    public void execute(final File file) {
                        SwingWorker sw = new SwingWorker() {
                            @Override
                            protected Object doInBackground() throws Exception {
                                HSSFWorkbook hssfWorkbook = tableModel.saveAsExcelSheet(new ProgressReporter() {
                                    @Override
                                    public void reportProgress(String note, int progress) {
                                        Vector v = new Vector();
                                        v.add(note);
                                        v.add(progress);
                                        publish(v);
                                    }
                                }, saveAsImage);
                                hssfWorkbook.write(new FileOutputStream(file));
                                hssfWorkbook.close();
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
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "Excel file saved.");
                                } catch (Exception e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                                }finally {
                                    progressMonitor.close();
                                }
                            }
                        };
                        sw.execute();

                    }
                });
            }
        });
        toolBar.add(excelBtn);

        JButton filterBtn = new JButton("Filter", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("filter-icon.png"))));
        filterBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                final FilterDialog dialog = new FilterDialog((JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor(),molecules.toArray(new PropertyMolecule[molecules.size()]));
                dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Boolean doInBackground() throws Exception {
                            PropertyMolecule[] selectedMolecules = dialog.getSelectedMolecules();
                            if(selectedMolecules!=null&&selectedMolecules.length>0) {
                                molecules.clear();
                                for (PropertyMolecule mol : selectedMolecules) {
                                    molecules.add(mol);
                                }
                                return new Boolean(true);
                            }
                            return new Boolean(false);
                        }

                        @Override
                        protected void done() {
                            try {
                                boolean result = (Boolean)get();
                                if(result) {
                                    updateMolTable();
                                }else{
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules selected.");
                                }
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                            }finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    sw.execute();
                }
            }
        });
        toolBar.add(filterBtn);
        JButton sdfBtn = new JButton("Export to 2D SDF", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("File-icon_small.png"))));
        sdfBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                    return;
                }
                FileUtil.saveToFile(currentDirectory, new FileNameExtensionFilter("SDF file","sdf"), new FileFunctor() {
                    @Override
                    public void execute(final File file) {
                        SwingWorker sw = new SwingWorker() {
                            @Override
                            protected Object doInBackground() throws Exception {
                                if(file!=null){
                                    File parent = file.getParentFile();
                                    if(parent!=null&&parent.exists()&&parent.isDirectory()){
                                        currentDirectory = parent;
                                    }
                                }
                                oemolostream ofs = new oemolostream();
                                ofs.SetFormat(OEFormat.SDF);
                                ofs.open(file.getAbsolutePath());
                                int progress = 0;
                                for(PropertyMolecule mol:molecules){
                                    progress ++;
                                    Vector v = new Vector();
                                    v.add(String.format("Saving molecule No. %d",progress));
                                    v.add(100*progress/molecules.size());
                                    publish(v);
                                    OEGraphMol mol1 = new OEGraphMol(mol.getMol());
                                    for(String propertyName:tableModel.selectedProperties){
                                        MolProperty property = mol.getProperty(propertyName);
                                        if(property!=null) {
                                            oechem.OESetSDData(mol1, propertyName, property.getProperty());
                                        }
                                    }
                                    oechem.OEWriteMolecule(ofs,mol1);
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
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "SDF file saved.");
                                } catch (Exception e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                                }finally {
                                    progressMonitor.close();
                                }
                            }
                        };
                        sw.execute();

                    }
                });
            }
        });
        toolBar.add(sdfBtn);

        JButton sdf3DBtn = new JButton("Export to 3D SDF", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("File-icon_small.png"))));
        sdf3DBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
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
                                for(PropertyMolecule mol:molecules){
                                    progress ++;
                                    Vector v = new Vector();
                                    v.add(String.format("Saving molecule No. %d",progress));
                                    v.add(100*progress/molecules.size());
                                    publish(v);
                                    OEGraphMol mol1 = new OEGraphMol(mol.getMol3d());
                                    for(String propertyName:tableModel.selectedProperties){
                                        MolProperty property = mol.getProperty(propertyName);
                                        if(property!=null) {
                                            oechem.OESetSDData(mol1, propertyName, property.getProperty());
                                        }
                                    }
                                    oechem.OEWriteMolecule(ofs,mol1);
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
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "SDF file saved.");
                                } catch (Exception e1) {
                                    e1.printStackTrace();
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                                }finally {
                                    progressMonitor.close();
                                }
                            }
                        };
                        sw.execute();

                    }
                });
            }
        });
        toolBar.add(sdf3DBtn);


        JButton vortexBtn = new JButton("Export to Vortex", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("vortex_icon.png"))));
        vortexBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                    return;
                }
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        File tempFile = File.createTempFile("Insilico", ".sdf");
                        oemolostream ofs = new oemolostream();
                        ofs.SetFormat(OEFormat.SDF);
                        ofs.open(tempFile.getAbsolutePath());
                        System.out.println(tempFile.getAbsolutePath());
                        int progress = 0;
                        for(PropertyMolecule mol:molecules){
                            progress ++;
                            Vector v = new Vector();
                            v.add(String.format("Saving molecule No. %d",progress));
                            v.add(100*progress/molecules.size());
                            publish(v);
                            OEGraphMol mol1 = new OEGraphMol(mol.getMol());
                            for(String propertyName:tableModel.selectedProperties){
                                MolProperty property = mol.getProperty(propertyName);
                                if(property!=null) {
                                    oechem.OESetSDData(mol1, propertyName, property.getProperty());
                                }
                            }
                            oechem.OEWriteMolecule(ofs,mol1);
                        }
                        ofs.close();
                        return tempFile;
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
                            File sdf = (File) get();
                            Desktop desktop = Desktop.getDesktop();
                            desktop.browse(new URI("http://10.74.2.128:8080/vortexweb/" + "vortex.jsp?file1=" + URLEncoder.encode(sdf.toURI().toURL().toString(), "UTF-8")));
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "File sent to Vortex.");
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                        }finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        });
        toolBar.add(vortexBtn);



        JButton calculateBtn = new JButton("Calculation", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("Calculator-icon.png"))));
        calculateBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                    return;
                }
                if(propertySelectionDialog==null){
                    propertySelectionDialog = new PropertySelectionDialog(ChemFunc.properties,ChemFunc.adme_models,ChemFunc.safety_models,ChemFunc.other_models);
                }
                propertySelectionDialog.setLocationRelativeTo(InSlilicoPanel.this);
                propertySelectionDialog.setVisible(true);
                if(propertySelectionDialog.isSubmitted){
                    final Vector<String> selectedProperties = propertySelectionDialog.getContentPane().getSelectedProperties();
                    String propertyString = Joiner.on(";").join(selectedProperties);
                    inSilicoTools.getInstance().logPrediction(propertyString,molecules.size());
                    String selectedTagTmp = null;
                    String unitTmp = null;
                    String selectedLogP = null;
                    if(selectedProperties.contains("Ligand Lipophilic Efficiency")){
                        LipophilicEfficiencyDialog dialog = new LipophilicEfficiencyDialog(existingTags);
                        dialog.setLocationRelativeTo(InSlilicoPanel.this);
                        dialog.setVisible(true);
                        if(dialog.isSubmitted()){
                            selectedLogP = dialog.getLogpType();
                            selectedTagTmp = dialog.getSelectedTag();
                            unitTmp = dialog.getUnit();
                        }
                    }
                    else if(selectedProperties.contains("Ligand Efficiency")){
                        LigandEfficiencyDialog dialog = new LigandEfficiencyDialog(existingTags);
                        dialog.setLocationRelativeTo(InSlilicoPanel.this);
                        dialog.setVisible(true);
                        if(dialog.isSubmitted()){
                            selectedTagTmp = dialog.getSelectedTag();
                            unitTmp = dialog.getUnit();
                        }
                    }
                    final String selectedTag = selectedTagTmp;
                    final String unit = unitTmp;
                    final String logpType = selectedLogP;
                    if(logpType!=null){
                        String newPropertyName = "Ligand Lipophilic Efficiency" + "_" + logpType;
                        selectedProperties.remove("Ligand Lipophilic Efficiency");
                        selectedProperties.add(newPropertyName);
                    }

                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            try {
//                                "MoKa LogP","MoKa LogD","MoKa pKa", "CNS MPO", "No. Aromatic Rings"
//                                "Passive Permeability","PgP Substrate", "Plasma Protein Binding (Human)",
//                                        "pH-Solubility Profile","pH-Lipophilicity Profile","VolSurf ADME Report","Structure Alerts","hERG inhibition"
//                                "cHLM","cHPPB","cMDR1","cRLM","cRPPB","cSolubility"
                                ChemFunc.calculateOEProperty(molecules);
                                if(selectedProperties.contains("cHLM")||selectedProperties.contains("cHPPB")||selectedProperties.contains("cMDR1")||selectedProperties.contains("cRLM")||selectedProperties.contains("cRPPB")||selectedProperties.contains("cSolubility")){
                                    ChemFunc.generateInSilicoADMEproperties(molecules);

                                }
                                if(selectedProperties.contains("Consensus LogP")||selectedProperties.contains("CNS PET MPO")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")||selectedProperties.contains("ChemAxon LogP")||(logpType!=null&&logpType.equals("ChemAxon LogP"))){
                                    ChemFunc.calculateChemAxonLogP(molecules,null);
                                }
                                if(selectedProperties.contains("Consensus LogP")||selectedProperties.contains("CNS mTEMPO")){
                                    for(PropertyMolecule mol:molecules){
                                        mol.addConsensusLogP();
                                    }
                                }
                                if(selectedProperties.contains("CNS mTEMPO")){
                                    ChemFunc.calculateCNSTempoScore(molecules);
                                }
                                if(selectedProperties.contains("BB Ratio")){
                                    ChemFunc.generateBBB(molecules);
                                }

                                if(selectedProperties.contains("ChemAxon LogD")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS PET MPO")){
                                    ChemFunc.calculateChemAxonLogDNeutral(molecules,null);
                                }
                                if(selectedProperties.contains("ChemAxon Acidic pKa")||selectedProperties.contains("ChemAxon Basic pKa")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS PET MPO")){
                                    ChemFunc.calculateChemAxon_pKa(molecules,null);
                                    if(!selectedProperties.contains("MoKa Basic pKa")&&!selectedProperties.contains("MoKa Acidic pKa")){
                                        InSilicoToolOptions.pka_type = PropertyMolecule.CHEMAXON_PKA;
                                    }
                                }

                                if(selectedProperties.contains("MoKa LogP") || selectedProperties.contains("MoKa LogD") || selectedProperties.contains("MoKa pKa")){
                                    ChemFunc.generateMoKaDescriptors(molecules);
                                }

                                if(selectedProperties.contains("CNS MPO")){
                                    ChemFunc.generateCNSMPODescriptors(molecules,false);
                                }

                                if(selectedProperties.contains("CNS PET MPO")){
                                    ChemFunc.generateCNSMPODescriptors(molecules,  true);
                                }


                                if(selectedProperties.contains("No. Aromatic Rings")){
                                    ChemFunc.calculateNoAroRings(molecules);
                                }

                                if(selectedProperties.contains("No. Aromatic Ring Systems")){
                                    ChemFunc.generateNoAromaticRingsDescriptors(molecules);
                                }

                                if(selectedProperties.contains("Plasma Protein Binding (Human)")){
                                    ChemFunc.generatePlasmaProteinBinding(molecules);
                                }
                                if(selectedProperties.contains("Structure Alerts")){
                                    StructureAlert[] alerts = ChemFunc.process_structure_alerts();
                                    for(PropertyMolecule mol:molecules){
                                        mol.addStructureAlerts(alerts);
                                    }
                                }
                                if(selectedProperties.contains("hERG inhibition")){
                                    ChemFunc.generateHERG(molecules);
                                }

                                if(selectedProperties.contains("RLM Qh%")){
                                    ChemFunc.generateRLM(molecules);
                                }

                                if(selectedProperties.contains("Efflux Ratio(B-A/A-B)1uM (SVM)")){
                                    ChemFunc.generateEffluxSVM(molecules);
                                }

                                if(selectedProperties.contains("Efflux Ratio(B-A/A-B)1uM (Kriging)")){
                                    ChemFunc.generateEffluxKrig(molecules);
                                }

                                if(selectedProperties.contains("Human VDss(L/kg)")){
                                    ChemFunc.generateVDss(molecules);
                                }

                                if(selectedProperties.contains("Ames")){
                                    ChemFunc.predictAMES(molecules);
                                }

                                if(selectedProperties.contains("Clearance Route")){
                                    ChemFunc.predictClearanceMechanism(molecules);
                                    selectedProperties.add("PC_1");
                                    selectedProperties.add("PC_2");
                                }

                                //"Metabolic Clearance","Renal Clearance"
                                if(selectedProperties.contains("Human Renal Clearance(mL/min/kg)")){
                                    ChemFunc.predictClearance(molecules,false);
                                }

                                if(selectedProperties.contains("Human Metabolic Clearance(mL/min/kg)")){
                                    ChemFunc.predictClearance(molecules,true);
                                }

                                if(selectedProperties.contains("Principal Moment of Inertia (NPR1,NPR2)")){
                                    ChemFunc.calculatePMIBatch(molecules);
                                    selectedProperties.remove("Principal Moment of Inertia (NPR1,NPR2)");
                                    selectedProperties.add("NPR1");
                                    selectedProperties.add("NPR2");
                                }

                                //"Volume3D", "SurfaceArea3D","PolarSurfaceArea3D", "DipoleMoment3D"
                                if(selectedProperties.contains("Volume3D")||selectedProperties.contains("SurfaceArea3D")
                                        ||selectedProperties.contains("PolarSurfaceArea3D")||selectedProperties.contains("DipoleMoment3D")){
                                    for(PropertyMolecule m:molecules) {
                                        OEChemFunc.getInstance().calculate3DSurfaceVolumeDipole(m);
                                    }
                                }

                                if(selectedProperties.contains("Metabolic Stability")){
                                    ChemFunc.calculateMetabolicStability(molecules);
                                }

                                if(selectedProperties.contains("Log_Solubility_pH7")){
                                    ChemFunc.calculateSolubility_pH_7(molecules);
                                }

                                if(selectedProperties.contains("HOMO")||selectedProperties.contains("LUMO")|selectedProperties.contains("HBS")){
                                    ChemFunc.getMopacProperties(molecules, new ProgressReporter() {
                                        @Override
                                        public void reportProgress(String note, int progress) {
                                            progressMonitor.setNote(note);
                                            progressMonitor.setProgress(progress);
                                        }
                                    });
                                }
                                if(selectedProperties.contains("hERG")){
                                    ChemFunc.predict_hERG_DL(molecules);
                                }

                                if(selectedProperties.contains("Papp A->B")){
                                    ChemFunc.predict_MDCK_DL(molecules);
                                }

                                if(selectedProperties.contains("General Metabolic Stability: T1/2 (min) (Mouse)")){
                                    ChemFunc.predict_MSM_DL(molecules);
                                }

                                if(selectedTag!=null&&unit!=null){
                                    double factor = 1e-6;
                                    if(unit.equals("nM")){
                                        factor = 1e-9;
                                    }
                                    for(PropertyMolecule mol:molecules){
                                        String tagData = oechem.OEGetSDData(mol.getMol(), selectedTag);
                                        System.out.println(tagData);
                                        if(tagData != null) {
                                            try {
                                                double value = Double.parseDouble(tagData);
                                                double pIC50 = - Math.log10(value*factor);
                                                if(logpType!=null) {
                                                    MolProperty property = mol.getProperty(logpType);
                                                    if(property!=null) {
                                                        double LLE = pIC50 - property.getValue();
                                                        String newPropertyName = "Ligand Lipophilic Efficiency" + "_" + logpType;
                                                        mol.addProperty(newPropertyName, "" + LLE);
                                                    }
                                                }
                                                double LE = 1.4 * pIC50 / mol.getHeavyAtomCount();
                                                mol.addProperty("Ligand Efficiency",""+LE);
                                            } catch (NumberFormatException e1) {

                                            }
                                        }
                                    }
                                }
                            } catch (IOException e1) {
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            }
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
                                for(PropertyMolecule p:molecules){
                                    OESDDataIter oesdDataIter = oechem.OEGetSDDataPairs(p.getMol());
                                    while(oesdDataIter.hasNext()){
                                        String tag = oesdDataIter.next().GetTag();
                                        if(!existingTags.contains(tag)){
                                            existingTags.add(tag);
                                        }
                                    }
                                }
                                for(String p:tableModel.getSelectedProperties()){
                                    if(existingTags.contains(p)){
                                        selectedProperties.add(p);
                                    }
                                }
                                tableModel.setSelectedProperties(selectedProperties);
                                updateMolTable();
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            }finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    sw.execute();
                }
            }
        });
        toolBar.add(calculateBtn);
        workspacePanel.add(toolBar, BorderLayout.NORTH);
    }

    private void updateStatusBar() {
        numMolLabel.setText(String.format("Number of Molecules (%d)", tableModel.propertyMolecules.size()));
    }

    private void buildTablePanel(){
        matrixMolTableModel = new MatrixMolTableModel(molecules);
        matrixTable = new MatrixCellTable(matrixMolTableModel);
        matrixTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(!e.getValueIsAdjusting()){
                    for(int i=0;i<matrixTable.getRowCount();i++){
                        for(int j=0;j<matrixTable.getColumnCount();j++){
                            Object value = matrixMolTableModel.getValueAt(i, j);
                            if(matrixTable.isCellSelected(i,j)){
                                if(value!=null&&value instanceof PropertyMolecule){
                                    ((PropertyMolecule)value).setIsSelected(true);
                                }
                            }else{
                                if(value!=null&&value instanceof PropertyMolecule){
                                    ((PropertyMolecule)value).setIsSelected(false);
                                }
                            }
                        }
                    }
                    tableModel.fireTableDataChanged();
                    matrixMolTableModel.fireTableDataChanged();
//                    for(Integer row:molTable.getSelectedRows()){
//                        int idx = molTable.convertRowIndexToModel(row);
//                        if(idx>=0){
//                            PropertyMolecule mol = molecules.get(idx);
//                            if(!mol.isSelected()){
//                                mol.setIsSelected(true);
//                            }
//                        }
//                    }
                }
            }
        });
//        matrixTable.setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer());
//        matrixTable.setShowGrid(true);
        matrixTable.setGridColor(Color.GRAY);
//        matrixTable.setRowHeight(150);
//        matrixTable.setShowGrid(true);
//        matrixTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        for(int i=0;i<matrixMolTableModel.getColumnCount();i++){
            matrixTable.getColumnModel().getColumn(i).setMinWidth(150);
        }

        tableModel = new PropertyMolTableModel(molecules,existingTags);
        molTable = new JXTable(tableModel);
        molTable.setShowGrid(true);
        molTable.setGridColor(Color.GRAY);
        molTable.setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer());
        molTable.setRowHeight(300);
        molTable.setSortOrderCycle(SortOrder.ASCENDING, SortOrder.DESCENDING, SortOrder.UNSORTED);
        molTable.setColumnControlVisible(true);
        molTable.setSearchable(new TableSearchable(molTable));
        molTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(!e.getValueIsAdjusting()){
                    for(Integer row:molTable.getSelectedRows()){
                        int idx = molTable.convertRowIndexToModel(row);
                        if(idx>=0){
                            PropertyMolecule mol = molecules.get(idx);
                            if(!mol.isSelected()){
                                mol.setIsSelected(true);
                            }
                        }
                    }
                    matrixMolTableModel.fireTableDataChanged();
                }
            }
        });
        molTable.getColumnModel().addColumnModelListener(new TableColumnModelListener() {
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
                if (molTable.getColumnModel().getColumnCount() > 0) {
                    int molColIdx = molTable.convertColumnIndexToModel(1);
                    molTableRowHeight = molTable.getColumnModel().getColumn(molColIdx).getWidth();
                    molTable.setRowHeight(molTableRowHeight);
                }
            }

            @Override
            public void columnSelectionChanged(ListSelectionEvent e) {

            }
        });
        buildPopupMenu();

        tabbedPane = new JideTabbedPane();
        tabbedPane.setShowCloseButton(true);
        add(tabbedPane, BorderLayout.CENTER);

        JPanel p1 = new JPanel(new BorderLayout());
        JToolBar bBar = new JToolBar();
        columnCountCB = new JComboBox(new Integer[]{2,3,4,5,6,7,8});
        columnCountCB.setSelectedItem(new Integer(4));
        columnCountCB.setMaximumSize(new Dimension(100,20));
        columnCountCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                matrixMolTableModel.setColumnCount((Integer)columnCountCB.getSelectedItem());
                fixMatrixTableFormat();
            }
        });
        bBar.add(new JToolBar.Separator());
        bBar.add(new JLabel("Column Count:"));
        bBar.add(columnCountCB);
        bBar.add(new JToolBar.Separator());
        bBar.add(new JLabel("Column Width:"));
        columnWidthSlider = new JSlider(0,250,0);
        columnWidthSlider.setMaximumSize(new Dimension(200,20));
        bBar.add(columnWidthSlider);
        columnWidthSlider.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                fixMatrixTableFormat();
            }
        });
        bBar.add(new JToolBar.Separator());
//        saveAsPdfBtn = new JButton("Save As Pdf");
//        saveAsPdfBtn.addActionListener(new ActionListener() {
//            @Override
//            public void actionPerformed(ActionEvent e) {
//                JFileChooser fileChooser = new JFileChooser();
//                FileNameExtensionFilter pdfFilter = new FileNameExtensionFilter("PDF (.pdf)", "pdf");
//                fileChooser.setFileFilter(pdfFilter);
//                if (fileChooser.showSaveDialog(InSlilicoPanel.this) == JFileChooser.APPROVE_OPTION) {
//                    try {
//                        File selectedFile = fileChooser.getSelectedFile();
//                        if(!selectedFile.getAbsolutePath().endsWith(".pdf")){
//                            selectedFile = new File(selectedFile+".pdf");
//                        }
//                        if(selectedFile.exists()){
//                            if(JOptionPane.showConfirmDialog(InSlilicoPanel.this,"File exists, overwrite?")==JOptionPane.NO_OPTION){
//                                return;
//                            }
//                        }
//                        int A4_width = 595;
//                        int A4_height = 842;
//                        int col = matrixMolTableModel.getColumnCount();
//                        int row = -1;
//                        if(col<=5){ //Portrait
//                            row = A4_height/(A4_width/col);
//                        }else{ //Landscape
//                            row = A4_width/(A4_height/col);
//                        }
//                        String result = ChemFunc.genMolPdfReport(selectedFile,molecules,row,col);
//                        if(result==null) {
//                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "PDF saved.");
//                        }else{
//                            JOptionPane.showMessageDialog(InSlilicoPanel.this, null);
//                        }
//                    } catch (Exception err) {
//                        err.printStackTrace();
//                        JOptionPane.showMessageDialog(InSlilicoPanel.this, err.getMessage());
//                    }
//                }
//            }
//        });
//        bBar.add(saveAsPdfBtn);

        p1.add(bBar,BorderLayout.SOUTH);
        p1.add(new JideScrollPane(matrixTable),BorderLayout.CENTER);

        workspaceTabbedPanel.add(new JideScrollPane(molTable),"Molecule List");
        workspaceTabbedPanel.add(p1,"Molecule Matrix");
        workspaceTabbedPanel.setTabClosableAt(0,false);
        workspaceTabbedPanel.setTabClosableAt(1,false);
    }

    class drawEganEggListener implements ActionListener{
        @Override
        public void actionPerformed(ActionEvent e) {
            if(molecules.size()==0){
                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                return;
            }
//            if(getSelectedMolecules().size()==0){
//                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules selected.");
//                return;
//            }
            StandardDialog dialog = new StandardDialog() {
                @Override
                public JComponent createBannerPanel() {
                    return null;
                }

                @Override
                public JComponent createContentPanel() {
                    String logPName = ChemFunc.getBestLogP(tableModel.selectedProperties);
                    if(plot==null){
                        plot = new EganEgg("<html><body><font color=\"#0000FF\">Absorption-CNS</font></body></html>", "2D tPSA", logPName, true, true, true);
                    }else{
                        plot.clear();
                        plot.getScatterplot().getPlot().getRangeAxis().setLabel(logPName);

                    }
                    int series = plot.addSeries("Molecule");
                    for(PropertyMolecule propertyMol:molecules){
                        try {
                            plot.addPoint(series, propertyMol.getProperty(PropertyMolecule.OE_TPSA).getValue(),propertyMol.getProperty(logPName).getValue(),propertyMol);
                        } catch (Exception e1) {
                            e1.printStackTrace();
                        }
                    }
                    return plot;
                }

                @Override
                public ButtonPanel createButtonPanel() {
                    ButtonPanel panel = new ButtonPanel();
                    panel.setAlignment(SwingConstants.CENTER);
                    JButton saveAsPdfBtn = new JButton("Save as Pdf",new ImageIcon(getClass().getClassLoader().getResource("PDF-icon_small.png")));
                    saveAsPdfBtn.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setDialogResult(RESULT_AFFIRMED);
                            setVisible(false);
                            if (plot != null) {
                                JFileChooser fileChooser = new JFileChooser();
                                FileNameExtensionFilter pdfFilter = new FileNameExtensionFilter("PDF (.pdf)", "pdf");
                                fileChooser.setFileFilter(pdfFilter);
                                if (fileChooser.showSaveDialog(InSlilicoPanel.this) == JFileChooser.APPROVE_OPTION) {
                                    try {
                                        File selectedFile = fileChooser.getSelectedFile();
                                        if(!selectedFile.getName().endsWith(".pdf")){
                                            selectedFile = new File(selectedFile+".pdf");
                                        }
                                        plot.writeChartAsPDF(new FileOutputStream(selectedFile));
                                        JOptionPane.showMessageDialog(InSlilicoPanel.this, "PDF saved.");
                                    } catch (Exception err) {
                                        err.printStackTrace();
                                        JOptionPane.showMessageDialog(InSlilicoPanel.this, err.getMessage());
                                    }
                                }
                            }
                            dispose();
                        }
                    });
                    panel.addButton(saveAsPdfBtn, ButtonPanel.AFFIRMATIVE_BUTTON);

                    JButton cancelBtn = new JButton("Cancel");
                    cancelBtn.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setDialogResult(RESULT_CANCELLED);
                            setVisible(false);
                            dispose();
                        }
                    });
                    panel.addButton(cancelBtn,ButtonPanel.CANCEL_BUTTON);

                    return panel;
                }
            };
            dialog.setModal(true);
            dialog.setAlwaysOnTop(true);
            dialog.setSize(new Dimension(800,600));
            dialog.setLocationRelativeTo(InSlilicoPanel.this);
            dialog.setVisible(true);
        }

    }

    class drawPMIListener implements ActionListener{
        @Override
        public void actionPerformed(ActionEvent e) {
            if(molecules.size()==0){
                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                return;
            }
//            if(getSelectedMolecules().size()==0){
//                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules selected.");
//                return;
//            }
            if((!tableModel.getSelectedProperties().contains("NPR1"))||(!tableModel.getSelectedProperties().contains("NPR2"))){
                JOptionPane.showMessageDialog(InSlilicoPanel.this,"Calculate Principal Moment of Inertia first.");
                return;
            }
            StandardDialog dialog = new StandardDialog() {
                @Override
                public JComponent createBannerPanel() {
                    return null;
                }

                @Override
                public JComponent createContentPanel() {
                    if(pmiPlot==null){
                        pmiPlot = new PrincipalMomentOfInertiaPlot("<html><body><font color=\"#0000FF\">Principal Moment of Inertia </font></body></html>", "NPR1", "NPR2", true, true);
                    }else{
                        pmiPlot.clear();

                    }
                    int series = pmiPlot.addSeries("Molecule");
                    for(PropertyMolecule propertyMol:molecules){
                        try {
                            MolProperty npr1 = propertyMol.getProperty("NPR1");
                            MolProperty npr2 = propertyMol.getProperty("NPR2");
                            if(npr1!=null&&npr2!=null) {
                                pmiPlot.addPoint(series, npr1.getValue(), npr2.getValue(), propertyMol);
                            }
                        } catch (Exception e1) {
                            e1.printStackTrace();
                        }
                    }
                    return pmiPlot;
                }

                @Override
                public ButtonPanel createButtonPanel() {
                    ButtonPanel panel = new ButtonPanel();
                    panel.setAlignment(SwingConstants.CENTER);
                    JButton saveAsPdfBtn = new JButton("Save as Pdf",new ImageIcon(getClass().getClassLoader().getResource("PDF-icon_small.png")));
                    saveAsPdfBtn.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setDialogResult(RESULT_AFFIRMED);
                            setVisible(false);
                            if (pmiPlot != null) {
                                JFileChooser fileChooser = new JFileChooser();
                                FileNameExtensionFilter pdfFilter = new FileNameExtensionFilter("PDF (.pdf)", "pdf");
                                fileChooser.setFileFilter(pdfFilter);
                                if (fileChooser.showSaveDialog(InSlilicoPanel.this) == JFileChooser.APPROVE_OPTION) {
                                    try {
                                        File selectedFile = fileChooser.getSelectedFile();
                                        if(!selectedFile.getName().endsWith(".pdf")){
                                            selectedFile = new File(selectedFile+".pdf");
                                        }
                                        pmiPlot.writeChartAsPDF(new FileOutputStream(selectedFile));
                                        JOptionPane.showMessageDialog(InSlilicoPanel.this, "PDF saved.");
                                    } catch (Exception err) {
                                        err.printStackTrace();
                                        JOptionPane.showMessageDialog(InSlilicoPanel.this, err.getMessage());
                                    }
                                }
                            }
                            dispose();
                        }
                    });
                    panel.addButton(saveAsPdfBtn, ButtonPanel.AFFIRMATIVE_BUTTON);

                    JButton cancelBtn = new JButton("Cancel");
                    cancelBtn.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setDialogResult(RESULT_CANCELLED);
                            setVisible(false);
                            dispose();
                        }
                    });
                    panel.addButton(cancelBtn,ButtonPanel.CANCEL_BUTTON);

                    return panel;
                }
            };
            dialog.setModal(true);
            dialog.setAlwaysOnTop(true);
            dialog.setSize(new Dimension(600,800));
            dialog.setLocationRelativeTo(InSlilicoPanel.this);
            dialog.setVisible(true);
        }

    }

    class drawClearanceMechanismListener implements ActionListener{
        @Override
        public void actionPerformed(ActionEvent e) {
            if(molecules.size()==0){
                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                return;
            }
//            if(getSelectedMolecules().size()==0){
//                JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules selected.");
//                return;
//            }
            if((!tableModel.getSelectedProperties().contains("Clearance Route"))){
                JOptionPane.showMessageDialog(InSlilicoPanel.this,"Run Prediction of Clearance Route first.");
                return;
            }
            StandardDialog dialog = new StandardDialog() {
                @Override
                public JComponent createBannerPanel() {
                    return null;
                }

                @Override
                public JComponent createContentPanel() {
                    if(clearancePlot==null){
                        clearancePlot = new ClearancePlot();
                    }else{
                        clearancePlot.clear();

                    }
                    int series = clearancePlot.addSeries("Molecule");
                    clearancePlot.getScatterplot().getPlot().getRenderer().setSeriesPaint(series,Color.BLACK);
                    for(PropertyMolecule propertyMol:molecules){
                        try {
                            MolProperty pc1 = propertyMol.getProperty("PC_1");
                            MolProperty pc2 = propertyMol.getProperty("PC_2");
                            if(pc1!=null&&pc2!=null) {
                                clearancePlot.addPoint(series, pc1.getValue(), pc2.getValue(), propertyMol);
                            }
                        } catch (Exception e1) {
                            e1.printStackTrace();
                        }
                    }
                    return clearancePlot;
                }

                @Override
                public ButtonPanel createButtonPanel() {
                    ButtonPanel panel = new ButtonPanel();
                    panel.setAlignment(SwingConstants.CENTER);
                    JButton saveAsPdfBtn = new JButton("Save as Pdf",new ImageIcon(getClass().getClassLoader().getResource("PDF-icon_small.png")));
                    saveAsPdfBtn.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setDialogResult(RESULT_AFFIRMED);
                            setVisible(false);
                            if (clearancePlot != null) {
                                JFileChooser fileChooser = new JFileChooser();
                                FileNameExtensionFilter pdfFilter = new FileNameExtensionFilter("PDF (.pdf)", "pdf");
                                fileChooser.setFileFilter(pdfFilter);
                                if (fileChooser.showSaveDialog(InSlilicoPanel.this) == JFileChooser.APPROVE_OPTION) {
                                    try {
                                        File selectedFile = fileChooser.getSelectedFile();
                                        if(!selectedFile.getName().endsWith(".pdf")){
                                            selectedFile = new File(selectedFile+".pdf");
                                        }
                                        clearancePlot.writeChartAsPDF(new FileOutputStream(selectedFile));
                                        JOptionPane.showMessageDialog(InSlilicoPanel.this, "PDF saved.");
                                    } catch (Exception err) {
                                        err.printStackTrace();
                                        JOptionPane.showMessageDialog(InSlilicoPanel.this, err.getMessage());
                                    }
                                }
                            }
                            dispose();
                        }
                    });
                    panel.addButton(saveAsPdfBtn, ButtonPanel.AFFIRMATIVE_BUTTON);

                    JButton cancelBtn = new JButton("Cancel");
                    cancelBtn.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setDialogResult(RESULT_CANCELLED);
                            setVisible(false);
                            dispose();
                        }
                    });
                    panel.addButton(cancelBtn,ButtonPanel.CANCEL_BUTTON);

                    return panel;
                }
            };
            dialog.setModal(true);
            dialog.setAlwaysOnTop(true);
            dialog.setSize(new Dimension(600,800));
            dialog.setLocationRelativeTo(InSlilicoPanel.this);
            dialog.setVisible(true);
        }

    }

    class LoadLibrarySessionListener implements ActionListener{
        public void actionPerformed(ActionEvent e){
            EnumerationWizard wizard = buildEnumerationWizard();

            JFileChooser fc = new JFileChooser(currentDirectory);
            fc.setMultiSelectionEnabled(false);
            fc.setFileFilter(new FileNameExtensionFilter("Library enumeration file.", "enum"));
            int option = fc.showOpenDialog(InSlilicoPanel.this);
            if (option == JFileChooser.APPROVE_OPTION) {
                currentDirectory = fc.getCurrentDirectory();
                final File file = fc.getSelectedFile();
                progressMonitor.setNote("Loading molecules ...");
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        if(!file.canRead()){
                            return null;
                        }
                        FileInputStream fin = new FileInputStream(file);
                        ObjectInputStream iip = new ObjectInputStream(fin);
                        String smirks = (String) iip.readObject();
                        ArrayList<ArrayList<ArrayList<SerializableMol>>> reagents = (ArrayList<ArrayList<ArrayList<SerializableMol>>>) iip.readObject();
                        ArrayList<SerializableMol> products = (ArrayList<SerializableMol>)iip.readObject();
                        Vector v = new Vector();
                        v.add(smirks);
                        v.add(reagents);
                        v.add(products);
                        iip.close();
                        fin.close();
                        return v;
                    }

                    @Override
                    protected void done() {
                        try {
                            Vector v = (Vector) get();
                            String smirks = (String) v.get(0);
                            ArrayList<ArrayList<ArrayList<SerializableMol>>> reagents = (ArrayList<ArrayList<ArrayList<SerializableMol>>>) v.get(1);
                            ArrayList<SerializableMol> products = (ArrayList<SerializableMol>)v.get(2);
                            if(smirks!=null&&smirks.length()>0&&reagents.size()>0&&products.size()>0){
                                wizard.loadSession(smirks,reagents,products);
                            }
                        } catch (Exception e1) {
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        }
    }

    class Load3DSessionListener implements ActionListener{
        public void actionPerformed(ActionEvent e){
            JFileChooser fc = new JFileChooser(currentDirectory);
            fc.setMultiSelectionEnabled(false);
            fc.setFileFilter(new FileNameExtensionFilter("3D Session file.", "ser"));
            int option = fc.showOpenDialog(InSlilicoPanel.this);
            if (option == JFileChooser.APPROVE_OPTION) {
                currentDirectory = fc.getCurrentDirectory();
                final File file = fc.getSelectedFile();
                progressMonitor.setNote("Loading molecules ...");
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        if(!file.canRead()){
                            return null;
                        }
                        FileInputStream fin = new FileInputStream(file);
                        ObjectInputStream iip = new ObjectInputStream(fin);
                        ArrayList<String> tags = (ArrayList<String>)iip.readObject();
                        ArrayList<SerializableMol> molList2 = (ArrayList<SerializableMol>) iip.readObject();
                        Vector v = new Vector();
                        v.add(tags);
                        v.add(molList2);
                        iip.close();
                        fin.close();
                        return v;
                    }

                    @Override
                    protected void done() {
                        try {
                            Vector v = (Vector) get();
                            ArrayList<String> tags = (ArrayList<String>)v.get(0);
                            ArrayList<SerializableMol> molList = (ArrayList<SerializableMol>) v.get(1);
                            if(molList!=null&&molList.size()>0){
                                Vector<PropertyMolecule> mols = new Vector<PropertyMolecule>();
                                String pdbStr = null;
                                String pdbName = null;
                                SerializableMol serializableMol = molList.get(0);
                                if (serializableMol.getType()==OEFormat.PDB){
                                    pdbStr = serializableMol.getMolStr();
                                    pdbName = serializableMol.getName();
                                    if(molList.size()>1) {
                                        for (int i = 1; i < molList.size(); i++) {
                                            OEGraphMol oeMol = molList.get(i).getOEMol();
                                            PropertyMolecule m = new PropertyMolecule(oeMol);
                                            for(String tag:tags) {
                                                if(oechem.OEHasSDData(oeMol,tag)) {
                                                    m.addProperty(tag, oechem.OEGetSDData(oeMol, tag));
                                                }
                                            }
                                            mols.add(m);
                                        }
                                    }
                                }else{
                                    for (int i = 0; i < molList.size(); i++) {
                                        OEGraphMol oeMol = molList.get(i).getOEMol();
                                        PropertyMolecule m = new PropertyMolecule(oeMol);
                                        for(String tag:tags) {
                                            if(oechem.OEHasSDData(oeMol,tag)) {
                                                m.addProperty(tag, oechem.OEGetSDData(oeMol, tag));
                                            }
                                        }
                                        mols.add(m);
                                    }
                                }
                                if(mols.size()>0){
                                    JFrame frame = new JFrame("Saved 3D Session");
                                    Mol3DTablePanel p = new Mol3DTablePanel(true);
                                    p.setPropertyMolecules(pdbName, pdbStr,mols);
                                    for(String tag:tags) {
                                        p.addProperty(tag);
                                    }
                                    frame.getContentPane().add(p);
                                    frame.setSize(new Dimension(1280,1024));
                                    frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                                    frame.setJMenuBar(p.getMenuBar(Mol3DTablePanel.DOCKING_MODE));
                                    frame.setVisible(true);
                                    RefineryUtilities.centerFrameOnScreen(frame);
                                    return;
                                }
                            }
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, "No molecule found.");
                        } catch (InterruptedException e1) {
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        }
    }




    class LoadSDFListener implements ActionListener{
        public void actionPerformed(ActionEvent e){
            JFileChooser fc = new JFileChooser(currentDirectory);
            fc.setMultiSelectionEnabled(true);
            fc.addChoosableFileFilter(new FileNameExtensionFilter("SDF file.", "sdf"));
            fc.addChoosableFileFilter(new FileNameExtensionFilter("ChemDraw CDX","cdx"));
            int option = fc.showOpenDialog(InSlilicoPanel.this);
            if (option == JFileChooser.APPROVE_OPTION) {
                currentDirectory = fc.getCurrentDirectory();
                final File[] files = fc.getSelectedFiles();
                molecules.clear();
                existingTags.clear();
//                    OEMolDatabase db = new OEMolDatabase(fname);
//                    final int numMols = db.GetMaxMolIdx();
                progressMonitor.setNote("Loading molecules ...");
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        for(File f :files) {
                            System.out.println(f.getAbsolutePath());
                            if(!f.canRead()){
                                continue;
                            }
                            oemolistream ifs = new oemolistream();
                            ifs.open(f.getAbsolutePath());
                            boolean isCdx = ifs.GetFormat()==OEFormat.CDX;
                            OEGraphMol mol = new OEGraphMol();
                            while (oechem.OEReadMolecule(ifs, mol)) {
                                if(isCdx){
                                    Vector<OEGraphMol> partMols = ChemFunc.splitMol(mol);
                                    for(OEGraphMol partMol:partMols){
                                        molecules.add(new PropertyMolecule(partMol));
                                    }
                                    break;
                                }
                                OESDDataIter oesdDataPairs = oechem.OEGetSDDataPairs(mol);
                                while (oesdDataPairs.hasNext()) {
                                    OESDDataPair p = oesdDataPairs.next();
                                    String tag = p.GetTag();
                                    if (!existingTags.contains(tag)) {
                                        existingTags.add(tag);
                                    }
                                }
                                molecules.add(new PropertyMolecule(mol));
                            }
                            ifs.close();
                        }
                        //ChemFunc.calculateOEProperty(molecules);
                        return molecules;
                    }

                    @Override
                    protected void done() {
                        try {
                            Vector<PropertyMolecule> molecules = (Vector<PropertyMolecule>) get();
                            molUpdateTable();
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, String.format("%d molecules loaded.", molecules.size()));
                        } catch (InterruptedException e1) {
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        }
    }

    protected void simSearchDb(String dbName) {
        JFrame frame = (JFrame) InSlilicoPanel.getInstance().getTopLevelAncestor();
        final SimilaritySearchDialog dialog = new SimilaritySearchDialog(frame,"Search "+dbName);
        dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
        dialog.setVisible(true);
        final Molecule mol = dialog.getMolecule();
        final float simCutoff = dialog.getSimCutoff();
        final int hits = dialog.getMaxNumHits();
        if(dialog.isCommitted()&& !mol.isEmpty()&&simCutoff>0&&simCutoff<=1.0){
            String smiles = null;
            try {
                smiles = MolExporter.exportToFormat(mol,"smiles:u,a");
            } catch (IOException e1) {
                e1.printStackTrace();
                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                return;
            }
            progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
            final String smiles1 = smiles;
            final SwingWorker checkWorker = new SwingWorker() {
                @Override
                protected Object doInBackground() throws Exception {
                    if(dbName.equalsIgnoreCase(MARKET_SELECT)) {
                        return FrontierDAO.getInstance().getSimilarMoleculeCountFromAldrich(smiles1, simCutoff, hits);
                    }else{
                        return InHouseCollectionDAO.getInstance().getSimilarMoleculeCountFromCellarity(smiles1,simCutoff,hits);
                    }
                }

                @Override
                protected void done() {
                    try {
                        int count = (Integer)get();
                        if(count ==0){
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"No Hits Found.");
                            progressMonitor.close();
                            return;
                        }
                        int option = JOptionPane.showConfirmDialog(InSlilicoPanel.this, String.format("%d hits found, continue?", count));
                        if(option!=JOptionPane.YES_OPTION){
                            progressMonitor.close();
                            return;
                        }
                        firePropertyChange("continue2search",true,false);
                    } catch (Exception e1) {
                        e1.printStackTrace();
                        progressMonitor.close();
                        JOptionPane.showMessageDialog(getInstance(),e1.getMessage());
                    }
                }

            };
            final SwingWorker searchWorker = new SwingWorker() {
                @Override
                protected Object doInBackground() throws Exception {
                    Vector<PropertyMolecule> propertyMolecules = null;
                    if(dbName.equalsIgnoreCase(MARKET_SELECT)){
                        propertyMolecules=FrontierDAO.getInstance().getSimilarMoleculesFromAldrich(smiles1, simCutoff, hits);
                    }else{
                        propertyMolecules=InHouseCollectionDAO.getInstance().getSimilarMoleculesFromCellarity(smiles1,simCutoff,hits);
                    }
                    ChemFunc.calculateOEProperty(propertyMolecules);
                    molecules.clear();
                    existingTags.clear();
                    molecules.addAll(propertyMolecules);
                    return propertyMolecules;
                }

                @Override
                protected void done() {
                    try {
                        Vector<PropertyMolecule> mols = (Vector<PropertyMolecule>)get();
                        molUpdateTable();
                        JOptionPane.showMessageDialog(InSlilicoPanel.this, String.format("%d molecules loaded.", mols.size()));
                    } catch (Exception e1) {
                        e1.printStackTrace();
                    }finally {
                        progressMonitor.close();
                    }
                }
            };
            checkWorker.addPropertyChangeListener(new PropertyChangeListener() {
                @Override
                public void propertyChange(PropertyChangeEvent evt) {
                    if(evt.getPropertyName().equals("continue2search")){
                        searchWorker.execute();
                    }
                }
            });
            checkWorker.execute();
        }
    }


    private void fixTableFormat(){
        if (tableModel.getColumnCount() > 0) {
            molTable.setRowHeight(molTableRowHeight);
            molTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
            int selectColIdx = molTable.convertColumnIndexToModel(0);
            int molColIdx = molTable.convertColumnIndexToModel(1);
            molTable.getColumnModel().getColumn(selectColIdx).setMaxWidth(20);
            molTable.getColumnModel().getColumn(molColIdx).setMinWidth(150);
        }
    }

    private void fixMatrixTableFormat(){
        matrixTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        for(int i=0;i<matrixMolTableModel.getColumnCount();i++){
            matrixTable.getColumnModel().getColumn(i).setPreferredWidth(150+columnWidthSlider.getValue());
        }
        matrixTable.setRowHeight(150+columnWidthSlider.getValue());
    }

    private void searchReagentDatabase(final String vendor){
        JFrame frame = (JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor();
        final SubstructureSearchDialog dialog ;
        if(vendor.equals(SKYHAWK)){
            dialog = new SubstructureSearchDialog(frame,"Search "+vendor, true, false);
        }else{
            dialog = new SubstructureSearchDialog(frame,"Search "+vendor);
        }
        dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
        dialog.setVisible(true);
        final String smarts = dialog.getSmarts();
        if(dialog.isCommitted()&&!smarts.isEmpty()){
            progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
            final SwingWorker checkWorker = new SwingWorker() {
                @Override
                protected Object doInBackground() throws Exception {
                    switch(vendor){
                        case ENAMINE:
                            int hits = dialog.getMaxNumHits();
                            float clogp = dialog.getCLogPLimit();
                            float mw = dialog.getMWLimit();
                            float psa = dialog.getPSALimit();
                            return FrontierDAO.getInstance().getCountsFromSmartsEnamine(smarts,hits, clogp, mw, psa);
                        case MARKET_SELECT:
                            hits = dialog.getMaxNumHits();
                            clogp = dialog.getCLogPLimit();
                            mw = dialog.getMWLimit();
                            psa = dialog.getPSALimit();
                            return FrontierDAO.getInstance().getCountsFromSmartsMarketSelect(smarts,hits,clogp,mw,psa);
                        case PHARMARON:
                            hits = dialog.getMaxNumHits();
                            clogp = dialog.getCLogPLimit();
                            mw = dialog.getMWLimit();
                            psa = dialog.getPSALimit();
                            return FrontierDAO.getInstance().getCountsFromSmartsPharmaron(smarts,hits,clogp,mw,psa);
                        case SKYHAWK:
                            hits = dialog.getMaxNumHits();
                            return InHouseCollectionDAO.getInstance().getCountFromSmartsFromCellarity(smarts,hits);
                        default:
                            throw new Exception("Unknown Vendor!");
                    }
                }

                @Override
                protected void done() {
                    try {
                        progressMonitor.close();
                        int count = (Integer)get();
                        if(count ==0){
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"No Hits Found.");
                            return;
                        }
                        int option = JOptionPane.showConfirmDialog(InSlilicoPanel.this, String.format("%d hits found, continue?", count));
                        if(option!=JOptionPane.YES_OPTION){
                            return;
                        }
                        firePropertyChange("continue2search",true,false);
                    } catch (Exception e1) {
                        e1.printStackTrace();
                    }
                }

            };
            final SwingWorker sw = new SwingWorker() {
                @Override
                protected Object doInBackground() throws Exception {
                    Vector<String> molFromSmarts;
                    switch(vendor){
                        case ENAMINE:
                            int hits = dialog.getMaxNumHits();
                            float clogp = dialog.getCLogPLimit();
                            float mw = dialog.getMWLimit();
                            float psa = dialog.getPSALimit();
                            molFromSmarts =  FrontierDAO.getInstance().getMolFromSmartsFromEnamine(smarts,hits, clogp, mw, psa);
                            break;
                        case MARKET_SELECT:
                            hits = dialog.getMaxNumHits();
                            clogp = dialog.getCLogPLimit();
                            mw = dialog.getMWLimit();
                            psa = dialog.getPSALimit();
                            molFromSmarts =  FrontierDAO.getInstance().getMolFromSmartsFromMarketSelect(smarts,hits, clogp, mw, psa);
                            break;
                        case PHARMARON:
                            hits = dialog.getMaxNumHits();
                            clogp = dialog.getCLogPLimit();
                            mw = dialog.getMWLimit();
                            psa = dialog.getPSALimit();
                            molFromSmarts = FrontierDAO.getInstance().getMolFromSmartsFromPharmaron(smarts,hits, clogp, mw, psa);
                            break;
                        case SKYHAWK:
                            hits = dialog.getMaxNumHits();
                            molFromSmarts = InHouseCollectionDAO.getInstance().getMolFromSmartsFromCellarity(smarts,hits);
                            break;
                        default:
                            throw new Exception("Unknown Vendor!");
                    }
                    Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
                    molecules.clear();
                    existingTags.clear();
                    for(String molString:molFromSmarts){
                        oemolistream ifs = new oemolistream();
                        ifs.SetFormat(OEFormat.SDF);
                        ifs.openstring(molString);
                        OEGraphMol oemol = new OEGraphMol();
                        oechem.OEReadMolecule(ifs, oemol);
                        ifs.close();
                        oechem.OESetSDData(oemol,"Name",oemol.GetTitle());
                        OESDDataIter oesdDataPairs = oechem.OEGetSDDataPairs(oemol);
                        while (oesdDataPairs.hasNext()) {
                            OESDDataPair p = oesdDataPairs.next();
                            String tag = p.GetTag();
                            if (!existingTags.contains(tag)) {
                                existingTags.add(tag);
                            }
                        }
                        propertyMolecules.add(new PropertyMolecule(oemol));
                    }
//                  ChemFunc.calculateOEProperty(propertyMolecules);
//                  ChemFunc.generateCLogP(propertyMolecules);
                    molecules.addAll(propertyMolecules);
                    return propertyMolecules;
                }

                @Override
                protected void done() {
                    try {
                        Vector<PropertyMolecule> mols = (Vector<PropertyMolecule>)get();
                        molUpdateTable();
                        progressMonitor.close();
                        JOptionPane.showMessageDialog(InSlilicoPanel.this, String.format("%d molecules loaded.", mols.size()));
                    } catch (Exception e1) {
                        e1.printStackTrace();
                    }finally {
                        progressMonitor.close();
                    }
                }
            };
            checkWorker.addPropertyChangeListener(new PropertyChangeListener() {
                @Override
                public void propertyChange(PropertyChangeEvent evt) {
                    if(evt.getPropertyName().equals("continue2search")){
                        progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                        sw.execute();
                    }
                }
            });
            checkWorker.execute();
        }



    }

    public JMenuBar getMenuBar(){
        JMenuBar menuBar = new JMenuBar();

        JMenu fileMenu = new JMenu("File");
        JMenuItem openSessionItem = new JMenuItem("Open Saved 3D Session");
        openSessionItem.addActionListener(new Load3DSessionListener());
        fileMenu.add(openSessionItem);

        JMenuItem openEnumSessionItem = new JMenuItem("Open Saved Library Enumeration Session");
        openEnumSessionItem.addActionListener(new LoadLibrarySessionListener());
        fileMenu.add(openEnumSessionItem);

        JMenuItem openItem = new JMenuItem("Open SDF file");
        openItem.addActionListener(new LoadSDFListener());
        fileMenu.add(openItem);

        final JMenuItem view3DItem = new JMenuItem("Open 3D Viewer");
        view3DItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFrame viewer3dframe = new JFrame("3D Viewer");
                Mol3DTablePanel panel = new Mol3DTablePanel();
                panel.setPreferredSize(new Dimension(1280,1024));
                if(!molecules.isEmpty()){
                    panel.setPropertyMolecules(molecules);
                }
                viewer3dframe.getContentPane().add(panel);
                viewer3dframe.setJMenuBar(panel.getMenuBar(Mol3DTablePanel.CONF_SEARCH_MODE));
                viewer3dframe.pack();
                viewer3dframe.setVisible(true);
            }
        });
        fileMenu.add(view3DItem);

        JMenuItem sketchItem = new JMenuItem("Sketch New Molecule");
        sketchItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                loadMolFromSketcher();
            }
        });
        fileMenu.add(sketchItem);

/*
        JMenuItem searchFrontierItem = new JMenuItem("Search Frontier Reagents ...");
        searchFrontierItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searchReagentDatabase(FRONTIER);
            }
        });
        fileMenu.add(searchFrontierItem);
*/
        JMenuItem searchEnamineItem = new JMenuItem("Search Enamine Reagents ...");
        searchEnamineItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searchReagentDatabase(ENAMINE);
            }
        });
        fileMenu.add(searchEnamineItem);

        JMenuItem searchPharmaronItem = new JMenuItem("Search Pharmaron Reagents ...");
        searchPharmaronItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searchReagentDatabase(PHARMARON);
            }
        });
        fileMenu.add(searchPharmaronItem);

        JMenuItem searchMarketSelectItem = new JMenuItem("Search Aldrich Market Select Reagents ...");
        searchMarketSelectItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searchReagentDatabase(MARKET_SELECT);
            }
        });
        fileMenu.add(searchMarketSelectItem);


        JMenuItem searchCellarityItem = new JMenuItem("Search Cellarity Compounds ...");
        searchCellarityItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searchReagentDatabase(SKYHAWK);
            }
        });
        fileMenu.add(searchCellarityItem);

        JMenuItem searchAldrichSimItem = new JMenuItem("Similarity Search Market Select ...");
        searchAldrichSimItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                simSearchDb(MARKET_SELECT);
            }
        });
        fileMenu.add(searchAldrichSimItem);

        JMenuItem searchCellaritySimItem = new JMenuItem("Similarity Search Cellarity ...");
        searchCellaritySimItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                simSearchDb(SKYHAWK);
            }
        });
        fileMenu.add(searchCellaritySimItem);



        /*
        JMenuItem searchEnaminePreferredItem = new JMenuItem("Search Enamine/My Reagents ...");
        searchEnaminePreferredItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                searchReagentDatabase(ENAMINE_PREFERRED);
            }
        });
        fileMenu.add(searchEnaminePreferredItem);

        JMenuItem searchGVKSimItem = new JMenuItem("Similarity Search GVK ...");
        searchGVKSimItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFrame frame = (JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor();
                final SimilaritySearchDialog dialog = new SimilaritySearchDialog(frame,"Search GVK");
                dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
                dialog.setVisible(true);
                final Molecule mol = dialog.getMolecule();
                final float simCutoff = dialog.getSimCutoff();
                final int hits = dialog.getMaxNumHits();
                if(dialog.isCommitted()&& !mol.isEmpty()&&simCutoff>0&&simCutoff<=1.0){
                    String smiles = null;
                    try {
                        smiles = MolExporter.exportToFormat(mol,"smiles:u,a");
                    } catch (IOException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                        return;
                    }
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    final String smiles1 = smiles;
                    final SwingWorker checkWorker = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            return FrontierDAO.getInstance().getSimilarMoleculeCountFromGVK(smiles1,simCutoff,hits);
                        }

                        @Override
                        protected void done() {
                            try {
                                int count = (Integer)get();
                                if(count ==0){
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No Hits Found.");
                                    progressMonitor.close();
                                    return;
                                }
                                int option = JOptionPane.showConfirmDialog(InSlilicoPanel.this, String.format("%d hits found, continue?", count));
                                if(option!=JOptionPane.YES_OPTION){
                                    progressMonitor.close();
                                    return;
                                }
                                firePropertyChange("continue2search",true,false);
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(getInstance(),e1.getMessage());
                            }finally {
                                progressMonitor.close();
                            }
                        }

                    };
                    final SwingWorker searchWorker = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            Vector<PropertyMolecule> propertyMolecules = FrontierDAO.getInstance().getSimilarMoleculesFromGVK(smiles1,simCutoff,hits);
                            ChemFunc.calculateOEProperty(propertyMolecules);
                            ChemFunc.generateCLogP(propertyMolecules);
                            molecules.clear();
                            existingTags.clear();
                            molecules.addAll(propertyMolecules);
                            return propertyMolecules;
                        }

                        @Override
                        protected void done() {
                            try {
                                Vector<PropertyMolecule> mols = (Vector<PropertyMolecule>)get();
                                molUpdateTable();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this, String.format("%d molecules loaded.", mols.size()));
                            } catch (Exception e1) {
                                e1.printStackTrace();
                            }finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    checkWorker.addPropertyChangeListener(new PropertyChangeListener() {
                        @Override
                        public void propertyChange(PropertyChangeEvent evt) {
                            if(evt.getPropertyName().equals("continue2search")){
                                searchWorker.execute();
                            }
                        }
                    });
                    checkWorker.execute();
                }

            }
        });
        fileMenu.add(searchGVKSimItem);

        JMenuItem searchChemblSimItem = new JMenuItem("Similarity Search ChemBl...");
        searchChemblSimItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFrame frame = (JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor();
                final SimilaritySearchDialog dialog = new SimilaritySearchDialog(frame,"Search ChemBl");
                dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
                dialog.setVisible(true);
                final Molecule mol = dialog.getMolecule();
                final float simCutoff = dialog.getSimCutoff();
                final int hits = dialog.getMaxNumHits();
                if(dialog.isCommitted()&& !mol.isEmpty()&&simCutoff>0&&simCutoff<=1.0){
                    String smiles = null;
                    try {
                        smiles = MolExporter.exportToFormat(mol,"smiles:u,a");
                    } catch (IOException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                        return;
                    }
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    final String smiles1 = smiles;
                    final SwingWorker checkWorker = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            return FrontierDAO.getInstance().getSimilarMoleculeCountFromChembl(smiles1,simCutoff,hits);
                        }

                        @Override
                        protected void done() {
                            try {
                                int count = (Integer)get();
                                if(count ==0){
                                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No Hits Found.");
                                    progressMonitor.close();
                                    return;
                                }
                                int option = JOptionPane.showConfirmDialog(InSlilicoPanel.this, String.format("%d hits found, continue?", count));
                                if(option!=JOptionPane.YES_OPTION){
                                    progressMonitor.close();
                                    return;
                                }
                                firePropertyChange("continue2search",true,false);
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                progressMonitor.close();
                                JOptionPane.showMessageDialog(getInstance(),e1.getMessage());
                            }
                        }

                    };
                    final SwingWorker searchWorker = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            Vector<PropertyMolecule> propertyMolecules = FrontierDAO.getInstance().getSimilarMoleculesFromChembl(smiles1,simCutoff,hits);
                            ChemFunc.calculateOEProperty(propertyMolecules);
                            ChemFunc.generateCLogP(propertyMolecules);
                            molecules.clear();
                            existingTags.clear();
                            molecules.addAll(propertyMolecules);
                            return propertyMolecules;
                        }

                        @Override
                        protected void done() {
                            try {
                                Vector<PropertyMolecule> mols = (Vector<PropertyMolecule>)get();
                                molUpdateTable();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this, String.format("%d molecules loaded.", mols.size()));
                            } catch (Exception e1) {
                                e1.printStackTrace();
                            }finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    checkWorker.addPropertyChangeListener(new PropertyChangeListener() {
                        @Override
                        public void propertyChange(PropertyChangeEvent evt) {
                            if(evt.getPropertyName().equals("continue2search")){
                                searchWorker.execute();
                            }
                        }
                    });
                    checkWorker.execute();
                }

            }
        });
        fileMenu.add(searchChemblSimItem);
        */


        JMenuItem enumerateItem = new JMenuItem("Enumeration from reaction");
        enumerateItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                inSilicoTools.getInstance().logEnumeration();
//                int tabIdx = -1;
//                for(int i=0;i<tabbedPane.getTabCount();i++){
//                    if(tabbedPane.getTitleAt(i).equals("Library Enumeration")){
//                        tabIdx = i;
//                        break;
//                    }
//                }
//                if(tabIdx >=0) {
//                    tabbedPane.setSelectedIndex(tabIdx);
//                    return;
//                }
                buildEnumerationWizard();
            }

        });
        fileMenu.add(enumerateItem);

        JMenuItem projectIdeaItem = new JMenuItem("View Project Idea");
        projectIdeaItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Vector<Project> projects = FrontierDAO.getInstance().getMy_projects();
                Project[] projectArray = projects.toArray(new Project[projects.size()]);
                Object selectedObj = JOptionPane.showInputDialog(InSlilicoPanel.this,
                        "Change project", "Project",
                        JOptionPane.QUESTION_MESSAGE,
                        null, projectArray,projectArray[0]);
                if(selectedObj!=null){
                    Project selectedProject = (Project) selectedObj;
                    String projectName = selectedProject.getProject_name();
                    inSilicoTools.getInstance().logCompoundTracking(projectName);
                    int tabIdx = -1;
                    for(int i=0;i<tabbedPane.getTabCount();i++){
                        if(tabbedPane.getTitleAt(i).equals("Project "+projectName)){
                            tabIdx = i;
                            break;
                        }
                    }
                    if(tabIdx>=0){
                        tabbedPane.setSelectedIndex(tabIdx);
                        return;
                    }
                    MyCompoundIdealPanel panel = new MyCompoundIdealPanel(selectedProject);
                    tabbedPane.addTab("Project "+projectName,panel);
                    tabbedPane.setSelectedIndex(tabbedPane.getTabCount()-1);
                }
            }
        });
        fileMenu.add(projectIdeaItem);

        JMenuItem coreItem = new JMenuItem("Register/View Cores");
        coreItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                try {
                    if(coresPanel==null){
                        Vector<Core> allMyCores = FrontierDAO.getInstance().getAllMyCores();
                        coresPanel = new MyCoreBrowsingPanel(allMyCores);
                    }
                    tabbedPane.addTab("Cores",coresPanel);
                    tabbedPane.setSelectedIndex(tabbedPane.getTabCount()-1);
                } catch (SQLException e) {
                    e.printStackTrace();
                }
            }
        });
        fileMenu.add(coreItem);


        JMenuItem purificationRequestItem = new JMenuItem("Purification request");
        purificationRequestItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(limsPanel==null){
                    limsPanel = new LIMSPanel();
                }
                tabbedPane.addTab("Purification Request", limsPanel);
                tabbedPane.setSelectedIndex(tabbedPane.getTabCount()-1);
            }
        });
//        fileMenu.add(purificationRequestItem);

        menuBar.add(fileMenu);

        JMenu editMenu = new JMenu("Edit");
        JMenuItem selectItem = new JMenuItem("Select all");
        selectItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (PropertyMolecule mol : molecules) {
                    mol.setIsSelected(true);
                }
                tableModel.fireTableDataChanged();
            }
        });
        editMenu.add(selectItem);


        JMenuItem unselectItem = new JMenuItem("Unselect all");
        unselectItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(PropertyMolecule mol:molecules){
                    mol.setIsSelected(false);
                }
                tableModel.fireTableDataChanged();
            }
        });
        editMenu.add(unselectItem);

        JMenuItem invertSelectionItem = new JMenuItem("Invert Selection");
        invertSelectionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (PropertyMolecule mol : molecules) {
                    if (mol.isSelected()) {
                        mol.setIsSelected(false);
                    } else {
                        mol.setIsSelected(true);
                    }
                }
                tableModel.fireTableDataChanged();
            }
        });
        editMenu.add(invertSelectionItem);

        JCheckBoxMenuItem showStructureAlertItem = new JCheckBoxMenuItem("Show Structure Alert",true);
        showStructureAlertItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JCheckBoxMenuItem b = (JCheckBoxMenuItem)e.getSource();
                if(b.isSelected()){
                    InSilicoToolOptions.show_structural_alert = true;
                }else{
                    InSilicoToolOptions.show_structural_alert = false;
                }
                updateMolTable();
            }
        });
        editMenu.add(showStructureAlertItem);

        JCheckBoxMenuItem useTrainedpKaItem = new JCheckBoxMenuItem("Use Trained ChemAxon pKa Model",false);
        useTrainedpKaItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JCheckBoxMenuItem b = (JCheckBoxMenuItem)e.getSource();
                if(b.isSelected()){
                    InSilicoToolOptions.use_trained_chemaxon_pKa = true;
                }else{
                    InSilicoToolOptions.use_trained_chemaxon_pKa = false;
                }
            }
        });
//        editMenu.add(useTrainedpKaItem);


        final JMenu pkaTypeMenu = new JMenu("pKa Type");
        final JMenu atomlabelTypeMenu = new JMenu("Show pKa/Charges");
        ButtonGroup bg1 = new ButtonGroup();
        JRadioButtonMenuItem showPartialChargesItem = new JRadioButtonMenuItem("Show charges on structure",false);
        showPartialChargesItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JRadioButtonMenuItem b = (JRadioButtonMenuItem)e.getSource();
                if(b.isSelected()){
                    InSilicoToolOptions.show_pKa = false;
                    InSilicoToolOptions.show_partial_charges=true;
                }else{
                    InSilicoToolOptions.show_pKa = true;
                    InSilicoToolOptions.show_pKa = false;
                }
                pkaTypeMenu.setEnabled(InSilicoToolOptions.show_pKa);
                updateMolTable();
            }
        });

        JRadioButtonMenuItem showPKaItem = new JRadioButtonMenuItem("Show pKa on structure",true);
        showPKaItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JRadioButtonMenuItem b = (JRadioButtonMenuItem)e.getSource();
                if(b.isSelected()){
                    InSilicoToolOptions.show_pKa = true;
                    InSilicoToolOptions.show_partial_charges=false;
                }else{
                    InSilicoToolOptions.show_pKa = false;
                    InSilicoToolOptions.show_pKa = true;
                }
                pkaTypeMenu.setEnabled(InSilicoToolOptions.show_pKa);
                updateMolTable();
            }
        });

        JRadioButtonMenuItem showNothingItem = new JRadioButtonMenuItem("No Label on structure", false);
        showNothingItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JRadioButtonMenuItem b = (JRadioButtonMenuItem)e.getSource();
                if(b.isSelected()) {
                    InSilicoToolOptions.show_pKa = false;
                    InSilicoToolOptions.show_partial_charges = false;
                }
                pkaTypeMenu.setEnabled(InSilicoToolOptions.show_pKa);
                updateMolTable();
            }
        });
        atomlabelTypeMenu.add(showNothingItem);
        atomlabelTypeMenu.add(showPartialChargesItem);
        atomlabelTypeMenu.add(showPKaItem);
        editMenu.add(atomlabelTypeMenu);
        bg1.add(showPartialChargesItem);
        bg1.add(showPKaItem);
        bg1.add(showNothingItem);

        ButtonGroup bg = new ButtonGroup();
        JRadioButtonMenuItem showPKaChemAxonItem = new JRadioButtonMenuItem("Show ChemAxon pKa",InSilicoToolOptions.pka_type==PropertyMolecule.CHEMAXON_PKA);
        showPKaChemAxonItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JRadioButtonMenuItem b = (JRadioButtonMenuItem) e.getSource();
                if(b.isSelected()){
                    InSilicoToolOptions.pka_type = PropertyMolecule.CHEMAXON_PKA;
                    if(InSilicoToolOptions.show_pKa) {
                        updateMolTable();
                    }
                }
            }
        });

        JRadioButtonMenuItem showPKaMoKaItem = new JRadioButtonMenuItem("Show MoKa pKa",InSilicoToolOptions.pka_type==PropertyMolecule.MOKA_PKA);
        showPKaMoKaItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JRadioButtonMenuItem b = (JRadioButtonMenuItem) e.getSource();
                if(b.isSelected()){
                    InSilicoToolOptions.pka_type = PropertyMolecule.MOKA_PKA;
                    if(InSilicoToolOptions.show_pKa) {
                        updateMolTable();
                    }
                }
            }
        });

        pkaTypeMenu.add(showPKaChemAxonItem);
        pkaTypeMenu.add(showPKaMoKaItem);

        editMenu.add(pkaTypeMenu);

        bg.add(showPKaChemAxonItem);
        bg.add(showPKaMoKaItem);

        JMenuItem clearAllItem = new JMenuItem("Clear All Molecules");
        clearAllItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int confirm = JOptionPane.showConfirmDialog(InSlilicoPanel.this,"All molecules are being removed, are you sure?");
                if(confirm==JOptionPane.YES_OPTION){
                    molecules.clear();
                    tableModel.selectedProperties.clear();
                    existingTags.clear();
                    tableModel.fireTableStructureChanged();
                    tableModel.fireTableDataChanged();
                    matrixMolTableModel.fireTableDataChanged();
                    updateStatusBar();
                }
            }
        });
        editMenu.add(clearAllItem);

        JMenuItem setNameItem = new JMenuItem("Set Molecule Name ...");
        setNameItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.isEmpty()||existingTags.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule or tag is available.");
                    return;
                }
                UserValuePickingDialog dialog = new UserValuePickingDialog(existingTags);
                dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
                dialog.setVisible(true);
                if(dialog.isSubmitted){
                    String tag = dialog.getSelectedTag();
                    for(PropertyMolecule m:molecules){
                        if(oechem.OEHasSDData(m.getMol(),tag)){
                            String title = oechem.OEGetSDData(m.getMol(), tag);
                            m.setName(title);
                        }else if (m.hasProperty(tag)) {
                            m.setName(m.getProperty(tag).getProperty());
                        }
                    }
                }
                updateMolTable();
            }
        });
        editMenu.add(setNameItem);

        JMenuItem setCASNOItem = new JMenuItem("Find CAS Number ...");
        setCASNOItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule is available.");
                    return;
                }
                if(!existingTags.contains("CAS_NO")){
                    existingTags.add("CAS_NO");
                }
                progressMonitor.setNote("Progress");
                progressMonitor.setMaximum(molecules.size());
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        FrontierDAO dao = FrontierDAO.getInstance();
                        int idx = 0;
                        for(PropertyMolecule m:molecules){
                            OEGraphMol mol = new OEGraphMol(m.getMol());
                            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                            String bio_number = null;
                            try {
                                bio_number = dao.getCASNumberFromStructure(oechem.OEMolToSmiles(mol));
                            } catch (SQLException e1) {
                                e1.printStackTrace();
                            }
                            if(bio_number!=null) {
                                m.addProperty("CAS_NO", bio_number);
                            }
                            publish(idx++);
                        }
                        return true;
                    }

                    @Override
                    protected void process(List chunks) {
                        int idx = (Integer)chunks.get(chunks.size() - 1);
                        progressMonitor.setProgress(idx);
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            Vector<String> selectedProperties = tableModel.getSelectedProperties();
                            if(!selectedProperties.contains("CAS_NO")) {
                                selectedProperties.add("CAS_NO");
                                tableModel.setSelectedProperties(selectedProperties);
                            }
                            updateMolTable();
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        });
        editMenu.add(setCASNOItem);


        JMenuItem findVIRNumberItem = new JMenuItem("Find VIR Number");
        findVIRNumberItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule is available.");
                    return;
                }
                if(!existingTags.contains("VIR-NUMBER")){
                    existingTags.add("VIR-NUMBER");
                }
                progressMonitor.setMillisToDecideToPopup(0);
                progressMonitor.setMillisToPopup(0);
                progressMonitor.setNote("Progress");
                progressMonitor.setMaximum(molecules.size());
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        InHouseCollectionDAO dao = InHouseCollectionDAO.getInstance();
                        int idx = 0;
                        for(PropertyMolecule m:molecules){
                            OEGraphMol mol = new OEGraphMol(m.getMol());
                            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                            String cy_number = null;
                            try {
                                String mol_string = ChemFunc.getCanonicalizedStructure(OEChemFunc.getInstance().getStringFromOEMol(mol));
                                cy_number = dao.getVIRNumberFromMol(mol_string);
                            } catch (SQLException e1) {
                                e1.printStackTrace();
                            }
                            if(cy_number!=null) {
                                m.addProperty("VIR-NUMBER", cy_number);
                            }
                            publish(idx++);
                        }
                        return true;
                    }

                    @Override
                    protected void process(List chunks) {
                        int idx = (Integer)chunks.get(chunks.size() - 1);
                        progressMonitor.setProgress(idx);
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            Vector<String> selectedProperties = tableModel.getSelectedProperties();
                            if(!selectedProperties.contains("VIR-NUMBER")) {
                                selectedProperties.add("VIR-NUMBER");
                                tableModel.setSelectedProperties(selectedProperties);
                            }
                            updateMolTable();
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        });
        editMenu.add(findVIRNumberItem);

        JMenuItem findCYNumberItem = new JMenuItem("Find CY Number ...");
        findCYNumberItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule is available.");
                    return;
                }
                if(!existingTags.contains("CY-NUMBER")){
                    existingTags.add("CY-NUMBER");
                }
                progressMonitor.setMillisToDecideToPopup(0);
                progressMonitor.setMillisToPopup(0);
                progressMonitor.setNote("Progress");
                progressMonitor.setMaximum(molecules.size());
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        InHouseCollectionDAO dao = InHouseCollectionDAO.getInstance();
                        int idx = 0;
                        for(PropertyMolecule m:molecules){
                            OEGraphMol mol = new OEGraphMol(m.getMol());
                            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                            String cy_number = null;
                            try {
                                String mol_string = ChemFunc.getCanonicalizedStructure(OEChemFunc.getInstance().getStringFromOEMol(mol));
                                cy_number = dao.getCYNumberFromMol(mol_string);
                            } catch (SQLException e1) {
                                e1.printStackTrace();
                            }
                            if(cy_number!=null) {
                                m.addProperty("CY-NUMBER", cy_number);
                            }
                            publish(idx++);
                        }
                        return true;
                    }

                    @Override
                    protected void process(List chunks) {
                        int idx = (Integer)chunks.get(chunks.size() - 1);
                        progressMonitor.setProgress(idx);
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            Vector<String> selectedProperties = tableModel.getSelectedProperties();
                            if(!selectedProperties.contains("CY-NUMBER")) {
                                selectedProperties.add("CY-NUMBER");
                                tableModel.setSelectedProperties(selectedProperties);
                            }
                            updateMolTable();
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        });
        editMenu.add(findCYNumberItem);

        JMenuItem setTautomerItem = new JMenuItem("Find CY-Number (Tautomer) ...");
        setTautomerItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule is available.");
                    return;
                }
                if(!existingTags.contains("CY-NUMBER(T)")){
                    existingTags.add("CY-NUMBER(T)");
                }
                progressMonitor.setMillisToDecideToPopup(0);
                progressMonitor.setMillisToPopup(0);
                progressMonitor.setNote("Progress");
                progressMonitor.setMaximum(molecules.size());
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        InHouseCollectionDAO dao = InHouseCollectionDAO.getInstance();
                        int idx = 0;
                        for(PropertyMolecule m:molecules){
                            OEGraphMol mol = new OEGraphMol(m.getMol());
                            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                            String bio_number = null;
                            try {
                                String mol_string = ChemFunc.getCanonicalizedStructure(OEChemFunc.getInstance().getStringFromOEMol(mol));
                                bio_number = dao.getCYTautomerFromMol(mol_string);
                            } catch (SQLException e1) {
                                e1.printStackTrace();
                            }
                            if(bio_number!=null) {
                                m.addProperty("CY-NUMBER(T)", bio_number);
                            }
                            publish(idx++);
                        }
                        return true;
                    }

                    @Override
                    protected void process(List chunks) {
                        int idx = (Integer)chunks.get(chunks.size() - 1);
                        progressMonitor.setProgress(idx);
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            Vector<String> selectedProperties = tableModel.getSelectedProperties();
                            if(!selectedProperties.contains("CY-NUMBER(T)")) {
                                selectedProperties.add("CY-NUMBER(T)");
                                tableModel.setSelectedProperties(selectedProperties);
                            }
                            updateMolTable();
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                        } finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        });
        editMenu.add(setTautomerItem);

        JMenuItem displayTagsItem = new JMenuItem("Toggle Tags ...");
        displayTagsItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.isEmpty()||existingTags.isEmpty()){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecule or tag is available.");
                    return;
                }
                String[] tags = existingTags.toArray(new String[existingTags.size()]);
                MultiChoicesDialog dialog = new MultiChoicesDialog((JFrame)InSlilicoPanel.this.getTopLevelAncestor(),tags, tableModel.getSelectedProperties());
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setVisible(true);
                if(dialog.isSubmitted()){
                    for(PropertyMolecule m:molecules){
                        for(String tag:dialog.getSelectedTags()){
                            if(oechem.OEHasSDData(m.getMol(), tag)){
                                m.addProperty(tag,oechem.OEGetSDData(m.getMol(),tag));
                            }
                        }
                    }
                    Vector<String> selectedProperties = tableModel.getSelectedProperties();
                    Vector<String> propertiesToRemove = new Vector<String>();
                    for(String tag:selectedProperties){
                        if(existingTags.contains(tag)){
                            propertiesToRemove.add(tag);
                        }
                    }
                    selectedProperties.removeAll(propertiesToRemove);
                    for(String tag:dialog.getSelectedTags()){
                        if(!selectedProperties.contains(tag)){
                            selectedProperties.add(tag);
                        }
                    }
                    tableModel.setSelectedProperties(selectedProperties);
                    updateMolTable();
                }
            }
        });
        editMenu.add(displayTagsItem);
        menuBar.add(editMenu);

        JMenu utilityMenu = new JMenu("Utilities");
        ImageIcon imageIcon = ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("Barcode-Scanner-icon.png")));
        JMenuItem barcodeScanItem = new JMenuItem("Find Tare Weight ...",imageIcon);
        barcodeScanItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                BarCodeScanDialog dialog = new BarCodeScanDialog();
                dialog.setSize(new Dimension(800,600));
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setVisible(true);
            }
        });
//        utilityMenu.add(barcodeScanItem);

        JMenuItem turnOnChiralMenu = new JMenuItem("Turn on MDL Chiral flags for selected molecules");
        turnOnChiralMenu.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(getSelectedMolecules().size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "No molecules are selected.");
                }
                ChemFunc.toggleMDLChiralFlags(molecules,true);
                updateMolTable();
            }
        });
        utilityMenu.add(turnOnChiralMenu);

        JMenuItem turnOffChiralMenu = new JMenuItem("Turn off MDL Chiral flags for selected molecules");
        turnOffChiralMenu.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(getSelectedMolecules().size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "No molecules are selected.");
                }
                ChemFunc.toggleMDLChiralFlags(molecules,false);
                updateMolTable();
            }
        });
        utilityMenu.add(turnOffChiralMenu);

        JMenuItem addOrChiralFlagMenuItem = new JMenuItem("Add OR flag for selected molecules ");
        addOrChiralFlagMenuItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(getSelectedMolecules().size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "No molecules are selected.");
                }
                ChemFunc.addOrChiralFlag(molecules);
                updateMolTable();
            }
        });
        utilityMenu.add(addOrChiralFlagMenuItem);

        JMenuItem deleteOrChiralFlagMenuItem = new JMenuItem("Delete OR flag for selected molecules ");
        deleteOrChiralFlagMenuItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(getSelectedMolecules().size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "No molecules are selected.");
                }
                ChemFunc.removeOrChiralFlags(molecules);
                updateMolTable();
            }
        });
        utilityMenu.add(deleteOrChiralFlagMenuItem);

        JMenuItem sdfComparatorItem = new JMenuItem("SDF Comparison/Filter");
        sdfComparatorItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                SDFComparisonTools tools = new SDFComparisonTools();
                JDialog dialog = new JDialog((JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor(), "SDF Comparison");
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.getContentPane().add(tools);
                dialog.pack();
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setResizable(false);
                dialog.setVisible(true);
            }
        });
        utilityMenu.add(sdfComparatorItem);

        JMenuItem nmrConversionItem = new JMenuItem("NMR pdb to lib ...");
        nmrConversionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                NMRConversionTool tools = new NMRConversionTool();
                JDialog dialog = new JDialog((JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor(), "NMR Conversion");
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.getContentPane().add(tools);
                dialog.pack();
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setResizable(false);
                dialog.setVisible(true);
            }
        });
        utilityMenu.add(nmrConversionItem);

        JMenuItem addTagItem = new JMenuItem("Add data tag...");
        addTagItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent actionEvent) {
                Vector<PropertyMolecule> selectedMolecules = getSelectedMolecules();
                if(selectedMolecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this, "No molecules are selected.");
                    return;
                }
                MultiInputTextDialog dialog = new MultiInputTextDialog((JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor());
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.setSize(new Dimension(300,180));
                dialog.setLocationRelativeTo(InSlilicoPanel.getInstance());
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    String tag = dialog.getTag();
                    String value = dialog.getValue();
                    for(PropertyMolecule mol: selectedMolecules){
                        oechem.OESetSDData(mol.getMol(),tag,value);
                    }
                }
            }
        });
        utilityMenu.add(addTagItem);

        JMenuItem chiralEnumerationItem = new JMenuItem("Chiral Enumeration");
        chiralEnumerationItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ChiralEnumerationTools tools = new ChiralEnumerationTools();
                JDialog dialog = new JDialog((JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor(), "Chiral Enumeration");
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.getContentPane().add(tools);
                dialog.pack();
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setResizable(false);
                dialog.setVisible(true);
            }
        });
//        utilityMenu.add(chiralEnumerationItem);


        JMenuItem predictKeyCmpdsMenu = new JMenuItem("Predict key compounds ...");
        predictKeyCmpdsMenu.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules!=null&&molecules.size()>0) {
                    if(molecules.size()>=3000){
                        JOptionPane.showMessageDialog(InSlilicoPanel.this,"Too many molecules (>=3000).");
                    }else {
                        HashMap<String,Number> defaultValues = new HashMap<String, Number>();
                        defaultValues.put("Cutoff",0.7);
                        defaultValues.put("Num. Hits", 5);
                        defaultValues.put("Max FP Radius",5);
                        HashMap<String,Range> rangeHashMap = new HashMap<String, Range>();
                        rangeHashMap.put("Cutoff",new NumericRange(0.5,1.0));
                        rangeHashMap.put("Num. Hits", new IntegerRange(1,10));
                        rangeHashMap.put("Max FP Radius", new IntegerRange(1,7));

                        MultiInputDialog dialog = new MultiInputDialog((JFrame)InSlilicoPanel.getInstance().getTopLevelAncestor(), defaultValues, rangeHashMap);
                        dialog.pack();
                        dialog.setLocationRelativeTo(InSlilicoPanel.this);
                        dialog.setVisible(true);
                        if(!dialog.isCommitted()){
                            return;
                        }
                        final HashMap<String, Number> resultMap = dialog.getResultMap();
                        SwingWorker sw = new SwingWorker() {
                            @Override
                            protected Object doInBackground() throws Exception {
                                ChemFunc.predictKeyCompounds(molecules, resultMap.get("Cutoff").floatValue(), resultMap.get("Num. Hits").intValue(), resultMap.get("Max FP Radius").intValue(), new ProgressReporter() {
                                    @Override
                                    public void reportProgress(final String note, final int progress) {
                                        SwingUtilities.invokeLater(new Runnable() {
                                            @Override
                                            public void run() {
                                                progressMonitor.setNote(note);
                                                progressMonitor.setProgress(progress);
                                            }
                                        });
                                    }
                                });
                                return null;
                            }

                            @Override
                            protected void done() {
                                try {
                                    get();
                                    if(!existingTags.contains("KeyCompound")) {
                                        existingTags.add("KeyCompound");
                                    }
                                    Vector<String> selectedProperties = tableModel.getSelectedProperties();
                                    if(!selectedProperties.contains("KeyCompound")) {
                                        selectedProperties.add("KeyCompound");
                                    }
                                    tableModel.setSelectedProperties(selectedProperties);
                                    updateMolTable();
                                } catch (InterruptedException e1) {
                                    e1.printStackTrace();
                                } catch (ExecutionException e1) {
                                    e1.printStackTrace();
                                } finally {
                                    progressMonitor.close();
                                }
                            }
                        };
                        sw.execute();
                    }
                }
            }
        });
        utilityMenu.add(predictKeyCmpdsMenu);

        JMenuItem volsurfItem = new JMenuItem("Generate VolSurf Descriptors ...");
        volsurfItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molecules.size()==0){
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,"No molecules available.");
                    return;
                }
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        ChemFunc.generateVolSurfDescriptors(molecules, tableModel.getSelectedProperties());
                        return null;
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            tableModel.setSelectedProperties(tableModel.getSelectedProperties());
                            for(String property:tableModel.getSelectedProperties()){
                                if(!existingTags.contains(property)){
                                    existingTags.add(property);
                                }
                            }
                            updateMolTable();
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"Process Interrupted.");
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(InSlilicoPanel.this,"Execution Exception.");
                        }finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();
            }
        });

//        utilityMenu.add(volsurfItem);

        JMenuItem checkStockItem = new JMenuItem("Evotec Inventory Check ...");
        checkStockItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                final String input = JOptionPane.showInputDialog("Lower Limit (\u00B5L):");
                if(input!=null){
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            double cutoff = Double.parseDouble(input);
                            String evotecInventoryInfo = ChemFunc.getEvotecInventoryInfo(cutoff);
                            return evotecInventoryInfo;
                        }

                        @Override
                        protected void done() {
                            try {
                                String resultStr = (String)get();
                                JSONParser parser = new JSONParser();
                                JSONObject jsonObject = (JSONObject) parser.parse(resultStr);
                                String csv = (String)jsonObject.get("csv");
                                long count = (Long)jsonObject.get("count");
                                File temp = File.createTempFile("tempfile", ".csv");
                                BufferedWriter bw = new BufferedWriter(new FileWriter(temp));
                                bw.write(csv);
                                bw.close();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,String.format("%d records found.",count));
                                Desktop.getDesktop().open(temp);
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            } catch (ParseException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            } catch (IOException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                            }finally {
                                progressMonitor.close();
                            }
                        }
                    };
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    sw.execute();
                }
                System.out.println(input);
            }
        });
//        utilityMenu.add(checkStockItem);

        menuBar.add(utilityMenu);
        JMenu modelingMenu = new JMenu("Modeling Application");
        JMenuItem confSearchItem = new JMenuItem("Conformational Search");
        confSearchItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int tabIdx = -1;
                for(int i=0;i<tabbedPane.getTabCount();i++){
                    if(tabbedPane.getTitleAt(i).equals("Conformation Search")){
                        tabIdx = i;
                        break;
                    }
                }
                if(tabIdx == -1) {
                    ConfSearchInputPanel panel = new ConfSearchInputPanel();
                    tabbedPane.add("Conformation Search", panel);
                    tabIdx = tabbedPane.getTabCount()-1;
                    tabbedPane.setTabClosableAt(tabIdx,true);
                    tabbedPane.setSelectedIndex(tabIdx);
                }else{
                    tabbedPane.setSelectedIndex(tabIdx);
                }
            }
        });
        modelingMenu.add(confSearchItem);

        JMenuItem superimposeItem = new JMenuItem("Superimpose");
        superimposeItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int tabIdx = -1;
                for(int i=0;i<tabbedPane.getTabCount();i++){
                    if(tabbedPane.getTitleAt(i).equals("Superimpose")){
                        tabIdx = i;
                        break;
                    }
                }
                if(tabIdx == -1) {
                    MolOverlayWizard panel = new MolOverlayWizard();
                    tabbedPane.add("Superimpose", panel);
                    tabIdx = tabbedPane.getTabCount()-1;
                    tabbedPane.setTabClosableAt(tabIdx,true);
                    tabbedPane.setSelectedIndex(tabIdx);
                }else{
                    tabbedPane.setSelectedIndex(tabIdx);
                }

            }
        });
        modelingMenu.add(superimposeItem);

        JMenuItem dockingItem = new JMenuItem("Docking");
        dockingItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int tabIdx = -1;
                for(int i=0;i<tabbedPane.getTabCount();i++){
                    if(tabbedPane.getTitleAt(i).equals("Docking")){
                        tabIdx = i;
                        break;
                    }
                }
                if(tabIdx == -1) {
                    MolDockingWizard panel = new MolDockingWizard();
                    tabbedPane.add("Docking", panel);
                    tabIdx = tabbedPane.getTabCount()-1;
                    tabbedPane.setTabClosableAt(tabIdx,true);
                    tabbedPane.setSelectedIndex(tabIdx);
                }else{
                    tabbedPane.setSelectedIndex(tabIdx);
                }
            }
        });
        modelingMenu.add(dockingItem);

        menuBar.add(modelingMenu);

        JMenu batchMenu = new JMenu("Batch jobs");
        JMenuItem descriptorCalculationItem = new JMenuItem("Calculate descriptors ...");
        descriptorCalculationItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JDialog dialog = new JDialog();
                BatchDescriptorGenerationPanel p = new BatchDescriptorGenerationPanel();
                dialog.getContentPane().add(p);
                dialog.setSize(new Dimension(720,800));
                dialog.setLocationRelativeTo(InSlilicoPanel.this);
                dialog.setVisible(true);
            }
        });
        batchMenu.add(descriptorCalculationItem);
        menuBar.add(batchMenu);

        JMenu plotMenu = new JMenu("Plot");
        JMenuItem eggItem = new JMenuItem("Egan Egg");
        eggItem.addActionListener(new drawEganEggListener());
        plotMenu.add(eggItem);

        JMenuItem pmiItem = new JMenuItem("Principal Moment of Inertia Plot");
        pmiItem.addActionListener(new drawPMIListener());
        plotMenu.add(pmiItem);

        JMenuItem clearanceRouteItem = new JMenuItem("Clearance Route Prediction Plot");
        clearanceRouteItem.addActionListener(new drawClearanceMechanismListener());
//        plotMenu.add(clearanceRouteItem);

        menuBar.add(plotMenu);
        JMenu helpMenu = new JMenu("Help"); //Should always be the last menuitem
        JMenuItem helpOnModel = new JMenuItem("Model Information");
        helpOnModel.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    Desktop desktop = Desktop.getDesktop();
                    desktop.browse(new URI("http://javelin.corp.My.com:8080/insilico/documentation/adme_properties.html"));
                } catch (Exception e1) {
                    JOptionPane.showMessageDialog(InSlilicoPanel.this,e1.getMessage());
                    e1.printStackTrace();
                }
            }
        });
        helpMenu.add(helpOnModel);
//        menuBar.add(helpMenu);

        return menuBar;
    }

    private void molUpdateTable() {
        tableModel.clearProperties();
        tableModel.fireTableDataChanged();
        tableModel.fireTableStructureChanged();
        fixTableFormat();
        updateStatusBar();

        matrixMolTableModel.fireTableDataChanged();
        matrixMolTableModel.fireTableStructureChanged();
        fixMatrixTableFormat();
    }

    private void updateMolTable() {
        if(molecules!=null){
            for(PropertyMolecule mol:molecules){
                mol.setNeedUpdate(true);
            }
        }
        tableModel.fireTableDataChanged();
        fixTableFormat();
        matrixMolTableModel.fireTableDataChanged();
        fixMatrixTableFormat();
        updateStatusBar();
    }

    private void loadMolFromSketcher() {
        sketcher.setLocationRelativeTo(this);
        sketcher.setVisible(true);
        if(sketcher.isCommitted()){
            try {
                Molecule molecule = sketcher.getMolecule();
                molecule.dearomatize();
                String molStr = MolExporter.exportToFormat(molecule,"sdf");
                OEGraphMol oemol = ChemFunc.getMolFromMolString(molStr, OEFormat.SDF);
                oemol.SetTitle(molecule.getName());
                if(oemol==null){
                    throw(new IOException("Failed to convert ChemAxon molecule to OEGraphMol."));
                }
                PropertyMolecule pmol = new PropertyMolecule(oemol);
                if(oemol.GetTitle()==null||oemol.GetTitle().isEmpty()){
                    oemol.SetTitle(pmol.getUniqName());
                }
                ChemFunc.calculateOEPropertySingleMol(pmol);
                molecules.add(pmol);
                if(molecules.size()==1){
                    tableModel.fireTableStructureChanged();
                }
                updateMolTable();
                updateStatusBar();
            } catch (IOException e1) {
                e1.printStackTrace();
                JOptionPane.showMessageDialog(InSlilicoPanel.getInstance(),e1.getMessage());
            }
        }
    }

    private EnumerationWizard buildEnumerationWizard() {
        int tabIdx;
        EnumerationWizard wizard = new EnumerationWizard();
        tabbedPane.add("Library Enumeration",wizard);
        tabIdx = tabbedPane.getTabCount()-1;
        tabbedPane.setTabClosableAt(tabIdx,true);
        tabbedPane.setSelectedIndex(tabIdx);
        wizard.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(evt.getPropertyName().equals("Finished")){
                    if(wizard.getProducts().isEmpty()){
                        return;
                    }
                    molecules.clear();
                    existingTags.clear();
                    progressMonitor.setNote("Loading molecules ...");
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            for(PropertyMolecule m:wizard.getProducts()){
                                OEGraphMol mol = m.getMol();
                                OESDDataIter oesdDataPairs = oechem.OEGetSDDataPairs(mol);
                                while(oesdDataPairs.hasNext()){
                                    OESDDataPair p = oesdDataPairs.next();
                                    String tag = p.GetTag();
                                    if(!existingTags.contains(tag)){
                                        existingTags.add(tag);
                                    }
                                }
                                molecules.add(new PropertyMolecule(mol));
                            }
                            SwingUtilities.invokeLater(new Runnable() {
                                @Override
                                public void run() {
                                    progressMonitor.setNote("Calculating descriptors...");
                                }
                            });
//                            ChemFunc.calculateOEProperty(molecules);
                            return molecules;
                        }

                        @Override
                        protected void done() {
                            try {
                                Vector<PropertyMolecule> molecules = (Vector<PropertyMolecule>) get();
                                molUpdateTable();
                                tabbedPane.setSelectedIndex(0);
                                JOptionPane.showMessageDialog(InSlilicoPanel.this, String.format("%d molecules loaded.", molecules.size()));
                            } catch (InterruptedException e1) {
                                JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                                e1.printStackTrace();
                            } catch (ExecutionException e1) {
                                JOptionPane.showMessageDialog(InSlilicoPanel.this, e1.getMessage());
                                e1.printStackTrace();
                            } finally {
                                progressMonitor.setNote("");
                                progressMonitor.close();
                            }
                        }
                    };
                    sw.execute();
                }
            }
        });
        return wizard;
    }


    public Vector<PropertyMolecule> getSelectedMolecules(){
        Vector<PropertyMolecule> selectedMolecules = new Vector<PropertyMolecule>();
        if(molecules!=null&&molecules.size()>0){
            for(PropertyMolecule molecule:molecules){
                if(molecule.isSelected()){
                    selectedMolecules.add(molecule);
                }
            }
        }
        return selectedMolecules;
    }

}
