package com.insilico.application.insilicotools.gui.widget;

import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.MolProperty;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.data.SerializableMol;
import com.insilico.application.insilicotools.gui.table.MolTable2D;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.gui.CompoundInputDialog;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.insilico.application.insilicotools.gui.PropertyMolTableModel;
import com.jidesoft.swing.JideScrollPane;
import com.jidesoft.swing.JideSplitPane;
import openeye.oechem.*;
import org.jfree.ui.RefineryUtilities;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 5/13/16.
 */
public class Mol3DTablePanel extends JPanel {
    Vector<PropertyMolecule> propertyMolecules;
    PropertyMolTableModel tableModel;
    Vector<String> propertyToDisplay;
    MolTable2D mol2dTable;
    MolViewer3D viewer3D;
    DesignProgressMonitor progressMonitor;
    CompoundInputDialog superimposeDialog;
    public static final int CONF_SEARCH_MODE = 1;
    public static final int SUPERIMPOSE_MODE = 2;
    public static final int DOCKING_MODE = 3;
    boolean useReceptor = false;
    String receptorPdbStr;
    String receptorName;
    HashMap<String,PropertyMolecule> molHash;

    public Mol3DTablePanel() {
        this(false);
    }


    public String getReceptorName() {
        return receptorName;
    }

    private void setReceptor(String receptorName, String receptor){
        this.receptorPdbStr = receptor;
        this.receptorName = receptorName;
        viewer3D.setReceptor(receptorName,receptor);
    }

    public void updateFromMolObject(MolObject molObject){
        if(molObject!=null&&molHash.containsKey(molObject.getName())){
            PropertyMolecule mol = molHash.get(molObject.getName());
            OEGraphMol mol3d = ChemFunc.getMolFromMolString(molObject.getMolString(), OEFormat.SDF);
            mol.setOEMol(mol3d);
            if(molObject.getGlide_score()!=null){
                mol.addProperty("Docking Score",molObject.getGlide_score().toString());
            }
        }
    }

    public Mol3DTablePanel(boolean useReceptor) {
        super(new BorderLayout());
        this.useReceptor = useReceptor;
        molHash = new HashMap<String, PropertyMolecule>();
        progressMonitor = new DesignProgressMonitor(Mol3DTablePanel.this,"","",0,100);
        viewer3D = new MolViewer3D();
        viewer3D.jymol.addPropertyChangeListener(new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                String molName = evt.getPropertyName();
                if(molHash.containsKey(molName)){
                    MolObject molObject = viewer3D.clickConsumer.molObjectHash.get(molName);
                    updateFromMolObject(molObject);
                }
            }
        });
        propertyMolecules = new Vector<PropertyMolecule>();
        propertyToDisplay = new Vector<String>();
        tableModel = new PropertyMolTableModel(propertyMolecules,propertyToDisplay,true);
        mol2dTable = new MolTable2D(tableModel);
        mol2dTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(!e.getValueIsAdjusting()){
                    for(PropertyMolecule mol:propertyMolecules){
                        if(mol.isSelected()){
                            viewer3D.addLigand(mol);
                        }else{
                            viewer3D.hideLigand(mol);
                        }
                    }
                    int[] selectedRows = mol2dTable.getSelectedRows();
                    for(Integer rowId:selectedRows){
                        PropertyMolecule mol = propertyMolecules.get(mol2dTable.convertRowIndexToModel(rowId));
                        viewer3D.addLigand(mol);
                    }
//                    tableModel.fireTableDataChanged();
                }
            }
        });


        JideSplitPane splitPane = new JideSplitPane(JideSplitPane.HORIZONTAL_SPLIT);
        splitPane.add(new JideScrollPane(mol2dTable));
        splitPane.add(viewer3D);
        splitPane.setShowGripper(true);
        splitPane.setProportionalLayout(true);
        splitPane.setProportions(new double[]{0.3});
        splitPane.setHeavyweightComponentEnabled(true);
        add(splitPane,BorderLayout.CENTER);
        add(viewer3D.buildBtnPanel(this.useReceptor),BorderLayout.SOUTH);

    }

    public void setPropertyMolecules(Vector<PropertyMolecule> propertyMolecules){
        viewer3D.setZoom(true);
        setPropertyMolecules(null, null,propertyMolecules);
    }

    public void setPropertyMolecules(String receptorName, String receptor, Vector<PropertyMolecule> propertyMolecules) {
        if(propertyMolecules==null||propertyMolecules.isEmpty()){
            this.propertyMolecules.clear();
            this.molHash.clear();
        }else{
            this.propertyMolecules.clear();
            this.propertyMolecules.addAll(propertyMolecules);
            for(PropertyMolecule mol:propertyMolecules){
                molHash.put(mol.getUniqName(),mol);
            }
        }
        mol2dTable.updateTable();
        viewer3D.clear();
        viewer3D.setZoom(true);
        if(receptor!=null) {
            setReceptor(receptorName, receptor);
        }
        if(this.propertyMolecules.size()>1) {
            this.propertyMolecules.get(0).setIsSelected(true);
            this.propertyMolecules.get(1).setIsSelected(true);
            mol2dTable.getSelectionModel().setSelectionInterval(0, 1);
        }else if(this.propertyMolecules.size()>0){
            this.propertyMolecules.get(0).setIsSelected(true);
            mol2dTable.getSelectionModel().setSelectionInterval(0, 0);
        }
    }

    public void addProperty(String property){
        Vector<String> selectedProperties = tableModel.getSelectedProperties();
        if(!selectedProperties.contains(property)) {
            selectedProperties.add(property);
            tableModel.setSelectedProperties(selectedProperties);
        }
    }

    public void removeProperty(String property){
        Vector<String> selectedProperties = tableModel.getSelectedProperties();
        if(selectedProperties.contains(property)) {
            selectedProperties.remove(property);
            tableModel.setSelectedProperties(selectedProperties);
        }
    }

    private void saveAsPymol(final boolean markedOnly) {
        int num  = 0;
        for(PropertyMolecule m:propertyMolecules){
            if(markedOnly){
                if(!m.isMarked()){
                    continue;
                }
            }else {
                if (!m.isSelected()) {
                    continue;
                }
            }
            num ++;
        }
        if(num == 0){
            JOptionPane.showMessageDialog(Mol3DTablePanel.this, "No molecules available.");
            return;
        }
        FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("Pymol Session file", "pse"), new FileFunctor() {
            @Override
            public void execute(final File file) {
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        JSONObject dict = new JSONObject();
                        if(useReceptor){
                            MolObject receptorObj = viewer3D.clickConsumer.getReceptorObj();
                            if(receptorObj!=null){
                                dict.put("receptor",receptorObj.getMolString());
                            }
                        }
                        JSONArray ligList = new JSONArray();
                        for(PropertyMolecule mol:propertyMolecules){
                            if(markedOnly){
                                if(!mol.isMarked()){
                                    continue;
                                }
                            }else{
                                if(!mol.isSelected()){
                                    continue;
                                }
                            }
                            String ligStr = ChemFunc.getMolString(mol.getMol3d());
                            if(ligStr!=null){
                                ligList.add(ligStr);
                            }
                        }

//                        Vector<MolObject> molObjects = viewer3D.clickConsumer.getMolObjects();
//                        if(molObjects!=null&&molObjects.size()>0) {
//                            for (MolObject molObject : molObjects) {
//                                if (molObject.isShowing()) {
//                                    ligList.add(molObject.getMolString());
//                                }
//                            }
//                        }
                        dict.put("ligands",ligList);
                        return (byte[]) ChemFunc.generatePymolSession(dict.toJSONString());
                    }

                    @Override
                    protected void done() {
                        progressMonitor.close();
                        try {
                            byte[] res = (byte[])get();
                            BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(file));
                            bos.write(res);
                            bos.close();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, "Pymol Session saved.");
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, e1.getMessage());
                        } finally {
                        }
                    }
                };
                sw.execute();
            }
        });
    }

    private void saveSession(final boolean markedOnly) {
        int num  = 0;
        for(PropertyMolecule m:propertyMolecules){
            num ++;
        }
        if(num == 0){
            JOptionPane.showMessageDialog(Mol3DTablePanel.this, "No molecules available.");
            return;
        }
        FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("Session file", "ser"), new FileFunctor() {
            @Override
            public void execute(final File file) {
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        ArrayList<String> tags = new ArrayList<String>();
                        for(String p:tableModel.getSelectedProperties()){
                            tags.add(p);
                        }

                        String outputfile = file.getAbsolutePath();
                        if(!outputfile.endsWith(".ser")){
                            outputfile = outputfile.split(".")[0]+".ser";
                        }
                        ArrayList<SerializableMol> molList = new ArrayList<SerializableMol>();
                        if(useReceptor){
                            SerializableMol receptorMol = new SerializableMol(receptorPdbStr, OEFormat.PDB);
                            receptorMol.setName(receptorName);
                            molList.add(receptorMol);
                        }
                        int progress = 0;
                        for (PropertyMolecule mol : propertyMolecules) {
                            if(markedOnly){
                                if(!mol.isMarked()){
                                    continue;
                                }
                            }
                            progress++;
                            Vector v = new Vector();
                            v.add(String.format("Saving molecule No. %d", progress));
                            v.add(100 * progress / propertyMolecules.size());
                            publish(v);
                            molList.add(mol.getSerializableMol());
                        }
                        FileOutputStream fos = new FileOutputStream(outputfile);
                        ObjectOutputStream oos = new ObjectOutputStream(fos);
                        oos.writeObject(tags);
                        oos.writeObject(molList);
                        oos.close();
                        fos.close();
                        return outputfile;
                    }

                    @Override
                    protected void process(java.util.List chunks) {
                        Vector v = (Vector) chunks.get(chunks.size() - 1);
                        String note = (String) v.get(0);
                        int progress = (Integer) v.get(1);
                        progressMonitor.setNote(note);
                        progressMonitor.setProgress(progress);
                    }

                    @Override
                    protected void done() {
                        progressMonitor.close();
                        try {
                            String output = (String)get();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, String.format("Session file saved as %s.",output));
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, e1.getMessage());
                        } finally {
                        }
                    }
                };
                sw.execute();

            }
        });
    }


    private void saveSDF(final boolean markedOnly) {
        int num  = 0;
        for(PropertyMolecule m:propertyMolecules){
            if(markedOnly){
                if(!m.isMarked()){
                    continue;
                }
            }else {
                if (!m.isSelected()) {
                    continue;
                }
            }
            num ++;
        }
        if(num == 0){
            JOptionPane.showMessageDialog(Mol3DTablePanel.this, "No molecules available.");
            return;
        }
        FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("SDF file", "sdf"), new FileFunctor() {
            @Override
            public void execute(final File file) {
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        oemolostream ofs = new oemolostream();
                        ofs.SetFormat(OEFormat.SDF);
                        ofs.open(file.getAbsolutePath());
                        int progress = 0;
                        for (PropertyMolecule mol : propertyMolecules) {
                            if(markedOnly){
                                if(!mol.isMarked()){
                                    continue;
                                }
                            }else {
                                if (!mol.isSelected()) {
                                    continue;
                                }
                            }
                            progress++;
                            Vector v = new Vector();
                            v.add(String.format("Saving molecule No. %d", progress));
                            v.add(100 * progress / propertyMolecules.size());
                            publish(v);
                            OEGraphMol mol1 = new OEGraphMol(mol.getMol3d());
                            for (String propertyName : tableModel.getSelectedProperties()) {
                                MolProperty property = mol.getProperty(propertyName);
                                if (property != null) {
                                    oechem.OESetSDData(mol1, propertyName, property.getProperty());
                                }
                            }
                            oechem.OEWriteMolecule(ofs, mol1);
                        }
                        ofs.close();
                        return null;
                    }

                    @Override
                    protected void process(java.util.List chunks) {
                        Vector v = (Vector) chunks.get(chunks.size() - 1);
                        String note = (String) v.get(0);
                        int progress = (Integer) v.get(1);
                        progressMonitor.setNote(note);
                        progressMonitor.setProgress(progress);
                    }

                    @Override
                    protected void done() {
                        progressMonitor.close();
                        try {
                            get();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, "SDF file saved.");
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, e1.getMessage());
                        } finally {
                        }
                    }
                };
                sw.execute();

            }
        });
    }

    public JMenuBar getMenuBar(int mode){
        JMenuBar menuBar = new JMenuBar();
        JMenu fileMenu = new JMenu("File");
        JMenuItem saveAsSessionItem = new JMenuItem("Save as session");
        saveAsSessionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveSession(false);
            }
        });
        fileMenu.add(saveAsSessionItem);

        JMenuItem saveMarkedAsSessionItem = new JMenuItem("Save marked as session");
        saveMarkedAsSessionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveSession(true);
            }
        });
        fileMenu.add(saveMarkedAsSessionItem);

        JMenuItem saveAsSdfMarkedItem = new JMenuItem("Save marked As Sdf");
        saveAsSdfMarkedItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveSDF(true);
            }

        });
        fileMenu.add(saveAsSdfMarkedItem);

        JMenuItem pymolMarkedItem = new JMenuItem("Save marked as pymol session ...");
        pymolMarkedItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveAsPymol(true);
            }
        });
        fileMenu.add(pymolMarkedItem);

        fileMenu.addSeparator();

        JMenuItem saveAsSdfItem = new JMenuItem("Save Selected As Sdf");
        saveAsSdfItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveSDF(false);
            }

        });
        fileMenu.add(saveAsSdfItem);

        JMenuItem pymolItem = new JMenuItem("Save selected as pymol session ...");
        pymolItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveAsPymol(false);
            }
        });
        fileMenu.add(pymolItem);

        JMenuItem loadPdbItem = new JMenuItem("Load Pdb ...");
        loadPdbItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser fc = new JFileChooser(InSlilicoPanel.getInstance().getCurrentDirectory());
                fc.setMultiSelectionEnabled(false);
                fc.setFileFilter(new FileNameExtensionFilter("Protein Data Bank file.", "pdb"));
                int option = fc.showOpenDialog(Mol3DTablePanel.this);
                if (option == JFileChooser.APPROVE_OPTION) {
                    InSlilicoPanel.getInstance().setCurrentDirectory(fc.getCurrentDirectory());
                    File file = fc.getSelectedFile();
                    if(file!=null&&file.canRead()){
                        try {
                            String pdbStr = FileUtil.readFileToString(file);
                            if(pdbStr!=null&&pdbStr.length()>0){
                                setReceptor("receptor",pdbStr);
                            }
                        } catch (IOException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(Mol3DTablePanel.this, e1.getMessage());
                        }
                    }else{
                        JOptionPane.showMessageDialog(Mol3DTablePanel.this, "Failed to open pdb file.");
                    }
                }
            }
        });
        fileMenu.add(loadPdbItem);

        JMenuItem freeformStrainEnergyItem = new JMenuItem("Calculate freeform strain energies (time consuming)");
        freeformStrainEnergyItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        ChemFunc.calculateFreeForm(propertyMolecules);
                        return null;
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            String[] propertyKeys = new String[]{"Erel","deltaG","Local Strain","Global Strain"};
                            for(String property:propertyKeys){
                                addProperty(property);
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
        });
        fileMenu.add(freeformStrainEnergyItem);



        menuBar.add(fileMenu);

        JMenu editMenu = new JMenu("Edit");
        JMenuItem pasteItem = new JMenuItem("Paste molecule 3D");
        pasteItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                try {
                    Object data = clipboard.getData(DataFlavor.stringFlavor);
                    if(data!=null){
                        String sdfString = (String)data;
                        OEGraphMol mol = ChemFunc.getMolFromMolString(sdfString, OEFormat.MDL);
                        if(mol!=null) {
                            if (mol.GetDimension() != 3) {
                                JOptionPane.showMessageDialog(Mol3DTablePanel.this,"Not a 3D molecule.");
                            }else{
                                PropertyMolecule propertyMolecule = new PropertyMolecule(mol);
                                propertyMolecules.add(propertyMolecule);
                                molHash.put(propertyMolecule.getUniqName(),propertyMolecule);
                                tableModel.fireTableDataChanged();
                            }
                        }

                    }
                } catch (UnsupportedFlavorException e1) {
                    e1.printStackTrace();
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
            }
        });
        editMenu.add(pasteItem);
        menuBar.add(editMenu);

        JMenu selectMenu = new JMenu("Select");

        JMenuItem markItem = new JMenuItem("Mark all");
        markItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (PropertyMolecule mol : propertyMolecules) {
                    mol.setMarked(true);
                }
                tableModel.fireTableDataChanged();
            }
        });
        selectMenu.add(markItem);


        JMenuItem unmarkItem = new JMenuItem("Unmark all");
        unmarkItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(PropertyMolecule mol:propertyMolecules){
                    mol.setMarked(false);
                }
                tableModel.fireTableDataChanged();
            }
        });
        selectMenu.add(unmarkItem);

        JMenuItem invertMarkItem = new JMenuItem("Invert Mark");
        invertMarkItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (PropertyMolecule mol : propertyMolecules) {
                    if (mol.isMarked()) {
                        mol.setMarked(false);
                    } else {
                        mol.setMarked(true);
                    }
                }
                tableModel.fireTableDataChanged();
            }
        });
        selectMenu.add(invertMarkItem);

        selectMenu.addSeparator();

        JMenuItem selectItem = new JMenuItem("Select all");
        selectItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (PropertyMolecule mol : propertyMolecules) {
                    mol.setIsSelected(true);
                    viewer3D.addLigand(mol);
                }
                tableModel.fireTableDataChanged();
            }
        });
        selectMenu.add(selectItem);


        JMenuItem unselectItem = new JMenuItem("Unselect all");
        unselectItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(PropertyMolecule mol:propertyMolecules){
                    mol.setIsSelected(false);
                    viewer3D.hideLigand(mol);
                }
                tableModel.fireTableDataChanged();
            }
        });
        selectMenu.add(unselectItem);

        JMenuItem invertSelectionItem = new JMenuItem("Invert Selection");
        invertSelectionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for (PropertyMolecule mol : propertyMolecules) {
                    if (mol.isSelected()) {
                        mol.setIsSelected(false);
                        viewer3D.hideLigand(mol);
                    } else {
                        mol.setIsSelected(true);
                        viewer3D.addLigand(mol);
                    }
                }
                tableModel.fireTableDataChanged();
            }
        });
        selectMenu.add(invertSelectionItem);

        JMenuItem selectHighlighedItem = new JMenuItem("Select highlighted");
        selectHighlighedItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(int idx:mol2dTable.getSelectedRows()){
                    if(idx>=0) {
                        int modelIdx = mol2dTable.convertRowIndexToModel(idx);
                        if(modelIdx>=0){
                            propertyMolecules.get(modelIdx).setIsSelected(true);
                        }
                    }
                }
            }
        });
        selectMenu.add(selectHighlighedItem);

        JMenuItem markHighlighedItem = new JMenuItem("Mark highlighted");
        markHighlighedItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(int idx:mol2dTable.getSelectedRows()){
                    if(idx>=0) {
                        int modelIdx = mol2dTable.convertRowIndexToModel(idx);
                        if(modelIdx>=0){
                            propertyMolecules.get(modelIdx).setMarked(true);
                        }
                    }
                }
            }
        });
        selectMenu.add(markHighlighedItem);




        menuBar.add(selectMenu);

        if(mode == CONF_SEARCH_MODE) {
            JMenu superimposeMenu = new JMenu("Superimpose");
            JMenuItem superimposeMenuItem = new JMenuItem("Superimpose");
            superimposeMenuItem.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    if (propertyMolecules == null || propertyMolecules.size() == 0) {
                        JOptionPane.showMessageDialog(Mol3DTablePanel.this, "No molecules available.");
                        return;
                    }

                    if (superimposeDialog == null) {
                        JFrame frame = (JFrame) Mol3DTablePanel.this.getTopLevelAncestor();
                        superimposeDialog = new CompoundInputDialog(frame);
                    }

                    superimposeDialog.setMolecule(propertyMolecules.get(0));
                    superimposeDialog.setLocationRelativeTo(Mol3DTablePanel.this);
                    superimposeDialog.setVisible(true);
                    if (superimposeDialog.isCommitted()) {
                        Molecule mol = superimposeDialog.getMolecule();
                        if (mol != null && !mol.isEmpty()) {
                            final OEGraphMol submol = OEChemFunc.getInstance().convertChemAxonMol(mol);
                            final OEGraphMol template = propertyMolecules.get(0).getMol3d();
                            if (template != null) {
                                SwingWorker sw = new SwingWorker() {
                                    @Override
                                    protected Object doInBackground() throws Exception {
                                        for (int i = 1; i < propertyMolecules.size(); i++) {
                                            publish(100*i/propertyMolecules.size());
                                            final PropertyMolecule pmol = propertyMolecules.get(i);
                                            OEGraphMol target3d = pmol.getMol3d();
                                            if (target3d == null) {
                                                continue;
                                            }
                                            Vector<Vector<Integer>> subSearchMatchingList = OEChemFunc.getInstance().getSubSearchMatchingList(submol, template, target3d);
                                            if (subSearchMatchingList != null && subSearchMatchingList.size() > 0) {
                                                OEGraphMol fittedMol = OEChemFunc.getInstance().matchByList(template, target3d, subSearchMatchingList);
                                                if (fittedMol != null) {
                                                    pmol.setOEMol3D(fittedMol);
                                                    SwingUtilities.invokeLater(new Runnable() {
                                                        @Override
                                                        public void run() {
                                                            viewer3D.updateLigand(pmol);
                                                        }
                                                    });
                                                }
                                            }

                                        }
                                        return null;
                                    }

                                    @Override
                                    protected void process(java.util.List chunks) {
                                        Integer progress = (Integer) chunks.get(chunks.size() - 1);
                                        progressMonitor.setNote("Align molecules ...");
                                        progressMonitor.setProgress(progress);
                                    }

                                    @Override
                                    protected void done() {
                                        try {
                                            get();
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
                    }
                }
            });
            superimposeMenu.add(superimposeMenuItem);
            menuBar.add(superimposeMenu);
        }

        return menuBar;
    }


    public static void main(String[] args) {
        JPopupMenu.setDefaultLightWeightPopupEnabled(false);
        try {
            File f = new File("/Users/jfeng1/oe_license.txt");
            oechem.OEAddLicenseData(FileUtil.readFileToString(f));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        oemolistream ifs = new oemolistream();
        ifs.open("/Users/jfeng1/all.sdf");
        OEGraphMol mol = new OEGraphMol();
        Vector<PropertyMolecule> molList = new Vector<PropertyMolecule>();
        while(oechem.OEReadMolecule(ifs,mol)){
            PropertyMolecule p = new PropertyMolecule(mol);
            p.addProperty("RMSD",oechem.OEGetSDData(mol,"RMSD"));
            p.addProperty("r_i_docking_score",oechem.OEGetSDData(mol,"r_i_docking_score"));
            molList.add(p);
        }
        ifs.close();

        ifs.open("/Users/jfeng1/00Demo/bart_receptor.pdb");
        OEGraphMol receptor = new OEGraphMol();
        oechem.OEReadMolecule(ifs,receptor);
        ifs.close();

        oemolostream ofs = new oemolostream();
        ofs.SetFormat(OEFormat.PDB);
        ofs.openstring();
        oechem.OEWriteMolecule(ofs,receptor);
        String pdbStr = ofs.GetString();

        JFrame frame = new JFrame("Test");
        Mol3DTablePanel p = new Mol3DTablePanel(true);
        p.setPropertyMolecules("bart_abeta",pdbStr,molList);
        p.addProperty("RMSD");
        p.addProperty("r_i_docking_score");
        frame.getContentPane().add(p);
        frame.setSize(new Dimension(1920,1800));
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setJMenuBar(p.getMenuBar(SUPERIMPOSE_MODE));
        frame.setVisible(true);
        RefineryUtilities.centerFrameOnScreen(frame);
    }
}
