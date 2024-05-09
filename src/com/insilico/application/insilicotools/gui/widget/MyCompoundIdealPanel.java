package com.insilico.application.insilicotools.gui.widget;

import chemaxon.license.LicenseManager;
import chemaxon.license.LicenseProcessingException;
import com.insilico.application.insilicotools.InSilicoToolOptions;
import com.insilico.application.insilicotools.data.*;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.insilico.application.insilicotools.gui.*;
import com.insilico.application.insilicotools.gui.table.MyCompoundTableModel;
import com.insilico.application.insilicotools.gui.table.MolTable2D;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.gui.util.TableUtils;
import com.insilico.application.insilicotools.inSilicoTools;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.util.OEChemWebLicenseInstaller;
import com.jgoodies.looks.plastic.Plastic3DLookAndFeel;
import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.plastic.theme.SkyBlue;
import com.jidesoft.swing.JideSplitPane;
import openeye.oechem.*;
import org.jdesktop.swingx.JXStatusBar;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.plaf.metal.MetalLookAndFeel;
import javax.swing.table.TableCellEditor;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.sql.Date;
import java.sql.SQLException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 7/28/17.
 */
public class MyCompoundIdealPanel extends JPanel {
    MyCompoundTableModel tableModel;
    MolTable2D table;
    Vector<Compound> MyCompounds;
    JButton addNewBtn;
    JButton insertBatchBtn ;
    JButton removeBtn;
    JButton editBtn;
    JButton commitBtn;
    DesignProgressMonitor progressMonitor = new DesignProgressMonitor(MyCompoundIdealPanel.this,"Progress","Progress",0,100);
    Project project;
    Vector<Compound> compoundsToRemove;
    JXStatusBar statusBar;
    JLabel numMolLabel;
    JTextPane textPane;
    JButton addCommentBtn;
    JButton deleteCommentBtn;
    Vector<Compound> originalCompounds;
    SubstructureSearchDialog substructureSearchDialog;
    boolean useSubStructureFilter = false;
    String substructure = null;
    Status statusDesired = null;
    HashMap<Integer,Boolean> substructureMap = new HashMap<>();
    JComboBox<String> operationCB = new JComboBox<>(new String[]{"Show","Hide"});
    JComboBox<String> chemistCB = new JComboBox();

    public MyCompoundIdealPanel(Project project) {
        super(new BorderLayout());
        this.compoundsToRemove = new Vector<>();
        this.project = project;
        this.originalCompounds = new Vector<>();
        this.MyCompounds = new Vector<>();
        statusBar = new JXStatusBar();
        statusBar.setPreferredSize(new Dimension(1000,40));
        numMolLabel = new JLabel("Number of Molecules (0), 0 deleted, 0 modified");
        statusBar.add(numMolLabel);
        tableModel = new MyCompoundTableModel(MyCompounds);
        tableModel.addTableModelListener(new TableModelListener() {
            @Override
            public void tableChanged(TableModelEvent e) {
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        int numModified = 0;
                        for(Compound c:MyCompounds){
                            if(c.isChanged()||c.isAdded()){
                                numModified += 1;
                            }
                        }
                        return numModified;
                    }

                    @Override
                    protected void done() {
                        try {
                            int numModified = (Integer)get();
                            numMolLabel.setText(String.format("Number of Molecules (%d),  %d deleted, %d modified.",MyCompounds.size(), compoundsToRemove.size(), numModified));
                        } catch (InterruptedException | ExecutionException e1) {
                            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                            e1.printStackTrace();
                        }
                    }
                };
                sw.execute();
            }
        });
        table = new MolTable2D(tableModel);
        int starIdx = table.convertColumnIndexToView(7);
        table.getColumn(starIdx).setCellRenderer(new StarRaterCellRender());
        table.getColumn(starIdx).setCellEditor(new StarRaterCellEditor());
        int viewIdx = table.convertColumnIndexToView(5);
        table.getColumn(viewIdx).setCellEditor(new StatusComboboxEditor(FrontierDAO.getInstance().getMy_status()));
        int assignedIdx = table.convertColumnIndexToModel(8);
        table.getColumn(assignedIdx).setCellEditor(new AssignedChemistComboboxEditor(FrontierDAO.getInstance().getChemistsByProject(project.getProject_id())));
        table.setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer());
        table.setRowHeight(100);
        table.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(table.getSelectedRowCount()!=1){
                    textPane.setText("");
                    addCommentBtn.setEnabled(false);
                    deleteCommentBtn.setEnabled(false);
                    return;
                }
                int row = table.getSelectedRow();
                if(row>=0&&table.convertRowIndexToModel(row)>=0){
                    addCommentBtn.setEnabled(true);
                    deleteCommentBtn.setEnabled(true);
                    int modelRow = table.convertRowIndexToModel(row);
                    Compound compound = MyCompounds.get(modelRow);
                    int mol_id = compound.getId();
                    try {
                        textPane.setText(FrontierDAO.getInstance().getMyCompoundComment(mol_id));
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                    }
                }
            }
        });
        JideSplitPane mainPanel = new JideSplitPane(JideSplitPane.VERTICAL_SPLIT);
        JPanel commentPanel = new JPanel(new BorderLayout());
        commentPanel.setBorder(new TitledBorder("Comment"));
        commentPanel.setPreferredSize(new Dimension(1280,100));
        textPane = new JTextPane();
        textPane.setEditable(false);
        JPanel commentBtnPanel = new JPanel();
        addCommentBtn = new JButton("Add Comment");
        String userName = InSlilicoPanel.getInstance().getUserName();
        addCommentBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(table.getSelectedRowCount()!=1){
                    //shouldn't happen.
                    return;
                }
                int row = table.getSelectedRow();
                if(row>=0&&table.convertRowIndexToModel(row)>=0){
                    int modelRow = table.convertRowIndexToModel(row);
                    Compound compound = MyCompounds.get(modelRow);
                    int mol_id = compound.getId();
                    String comment = JOptionPane.showInputDialog("");
                    if(comment.trim().isEmpty()){
                        return;
                    }
                    try {
                        int chemist_id = FrontierDAO.getInstance().getChemistId(userName);
                        FrontierDAO.getInstance().addComment(mol_id,chemist_id,comment);
                        textPane.setText("");
                        textPane.setText(FrontierDAO.getInstance().getMyCompoundComment(mol_id));
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                    }
                }
            }
        });
        deleteCommentBtn = new JButton("Delete my comments.");
        deleteCommentBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(table.getSelectedRowCount()!=1){
                    //shouldn't happen.
                    return;
                }
                int row = table.getSelectedRow();
                if(row>=0&&table.convertRowIndexToModel(row)>=0){
                    int modelRow = table.convertRowIndexToModel(row);
                    Compound compound = MyCompounds.get(modelRow);
                    int mol_id = compound.getId();
                    try {
                        int chemist_id = FrontierDAO.getInstance().getChemistId(userName);
                        FrontierDAO.getInstance().deleteComment(mol_id, chemist_id);
                        textPane.setText("");
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                    }
                }
            }
        });
        commentBtnPanel.add(addCommentBtn);
        commentBtnPanel.add(deleteCommentBtn);
        commentPanel.add(new JScrollPane(textPane),BorderLayout.CENTER);
        commentPanel.add(commentBtnPanel,BorderLayout.SOUTH);
        mainPanel.addPane(new JScrollPane(table));
        mainPanel.addPane(commentPanel);
        mainPanel.setProportionalLayout(true);
        mainPanel.setProportions(new double[]{0.8});
        add(mainPanel,BorderLayout.CENTER);
        add(buildBtnPanel(),BorderLayout.NORTH);
        add(statusBar,BorderLayout.SOUTH);

        try {
            setProject(project);
        } catch (SQLException e) {
            e.printStackTrace();
            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e.getMessage());
        }
    }

    private void calculateProperties(final List<Compound> cmpdList) throws Exception {
        final Vector<String> selectedProperties = new Vector<>();
        selectedProperties.addAll(Arrays.asList(Compound.PROPERTIES_TO_USE));
        final Vector<PropertyMolecule> molecules = new Vector<>();
        for(Compound c:cmpdList){
            if(c.hasProperty()){
                continue;
            }
            molecules.add(c.getPropertyMol());
        }
        if(molecules.isEmpty()){
            return;
        }
        ChemFunc.calculateOEProperty(molecules);
        if(selectedProperties.contains("Consensus LogP")||selectedProperties.contains("CNS mTEMPO")||selectedProperties.contains("ChemAxon LogP")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")){
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

        if(selectedProperties.contains("ChemAxon LogD")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")){
            ChemFunc.calculateChemAxonLogDNeutral(molecules,null);
        }
        if(selectedProperties.contains("ChemAxon Acidic pKa")||selectedProperties.contains("ChemAxon Basic pKa")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")){
            ChemFunc.calculateChemAxon_pKa(molecules,null);
            if(!selectedProperties.contains("MoKa Basic pKa")&&!selectedProperties.contains("MoKa Acidic pKa")){
                InSilicoToolOptions.pka_type = PropertyMolecule.CHEMAXON_PKA;
            }
        }
        if(selectedProperties.contains("CNS MPO")){
            ChemFunc.generateCNSMPODescriptors(molecules, false);
        }

        if(selectedProperties.contains("No. Aromatic Rings")){
            ChemFunc.calculateNoAroRings(molecules);
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

        if(selectedProperties.contains("Efflux Ratio(B-A/A-B)1uM")){
            ChemFunc.generateEffluxKrig(molecules);
            ChemFunc.generateEffluxSVM(molecules);
        }


        if(selectedProperties.contains("VDss(L/kg)")){
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
        for(Compound c:cmpdList){
            if(!c.hasProperty()){
                c.save();
            }
        }

    }

    private void findBioNumbers(final List<Compound> cmpdList){
        SwingWorker sw = new SwingWorker() {
            @Override
            protected Object doInBackground() throws Exception {
                InHouseCollectionDAO dao = InHouseCollectionDAO.getInstance();
                final Vector<PropertyMolecule> molecules = new Vector<>();
                for(Compound c:cmpdList){
                    molecules.add(c.getPropertyMol());
                }

                for(PropertyMolecule m:molecules){
                    OEGraphMol mol = new OEGraphMol(m.getMol());
                    oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                    String cy_number = null;
                    String cy_number_t = null;
                    try {
                        cy_number = dao.getCYNumberFromMol(OEChemFunc.getInstance().getStringFromOEMol(mol));
                        cy_number_t = dao.getCYTautomerFromMol(OEChemFunc.getInstance().getStringFromOEMol(mol));
                    } catch (SQLException e1) {
                        e1.printStackTrace();
                    }
                    if(cy_number!=null) {
                        m.addProperty("CY-Number", cy_number);
                    }
                    if(cy_number_t!=null) {
                        m.addProperty("CY-Number(T)", cy_number_t);
                    }
                }
                return null;
            }

            @Override
            protected void done() {
                try {
                    get();
                    tableModel.fireTableDataChanged();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    e.printStackTrace();
                }
            }
        };
        sw.execute();
    }

    private JToolBar buildBtnPanel(){
        JToolBar p = new JToolBar();
        addNewBtn = new JButton("Add New Idea ...");
        addNewBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MyCompoundModifyDialog dialog = new MyCompoundModifyDialog(project);
                dialog.setLocationRelativeTo(MyCompoundIdealPanel.this);
                dialog.setVisible(true);
                Compound MyCompound = dialog.getMyCompound();
                if(MyCompound !=null){
                    try {
                        calculateProperties(Collections.singletonList(MyCompound));
                        findBioNumbers(Collections.singletonList(MyCompound));
                        int id = MyCompound.save();
                        if(id==-1){
                            throw new Exception("Failed to save compound, duplicate?");
                        }
                        MyCompounds.add(MyCompound);
                        originalCompounds.add(MyCompound);
                        tableModel.fireTableRowsInserted(MyCompounds.size()-1,MyCompounds.size()-1);
                    } catch (Exception e1) {
                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                        e1.printStackTrace();
                        return;
                    }
                }
            }
        });

        insertBatchBtn = new JButton("Batch Insert Ideas ...");
        String userName = InSlilicoPanel.getInstance().getUserName();
        insertBatchBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                final MyCompoundsBatchRegistrationDialog dialog = new MyCompoundsBatchRegistrationDialog();
                dialog.setMyProject(project);
                dialog.setLocationRelativeTo(MyCompoundIdealPanel.this);
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    final String sdfPathName = dialog.getSdfFilePath();
                    final Project project = dialog.getMyProject();
                    final Status status = dialog.getMyStatus();
                    final String nameTag = dialog.getMolNameTag();
                    final int chemist_id = FrontierDAO.getInstance().getChemistId(userName);
                    if(status==null){
                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,"No status specified");
                        return;
                    }
                    if(project==null){
                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,"No project specified");
                        return;
                    }

                    progressMonitor.setMillisToDecideToPopup(0);
                    progressMonitor.setMillisToPopup(0);
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            Vector<Integer> idList = FrontierDAO.getInstance().insertMyCompoundsBatch(sdfPathName, nameTag, project.getProject_id(), status.getStatus_id(), chemist_id);
                            return idList;
                        }

                        @Override
                        protected void done() {
                            try {
                                Vector<Integer> idList = (Vector<Integer>)get();
                                if(idList.size()>0) {
                                    reload();
                                }
                                progressMonitor.close();
                                JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,String.format("%d molecules imported.",idList.size()));
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
        removeBtn = new JButton("Remove Selected Ideas ...");
        removeBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int tryToRemoveOthersIdea = 0;
                for(int i=MyCompounds.size()-1;i>-1;i--){
                    Compound c = MyCompounds.get(i);
                    if(c.getSelected()){
                        if(!c.getChemist().equals("anonymous") &&!c.getChemist().equals(userName)){
                            tryToRemoveOthersIdea+=1;
                        }else {
                            compoundsToRemove.add(c);
                        }
                    }
                }
                MyCompounds.removeAll(compoundsToRemove);
                originalCompounds.removeAll(compoundsToRemove);
                tableModel.fireTableDataChanged();
                if(tryToRemoveOthersIdea>0) {
                    JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,String.format("%d ideas are from other people and will not be removed.",tryToRemoveOthersIdea));
                }
            }
        });
        editBtn = new JButton("Change Current Idea ...");
        editBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int row = table.getSelectedRow();
                if(row<0){
                    JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,"No molecule selected!");
                    return;
                }
                int row_model = table.convertRowIndexToModel(row);
                if(row_model>=0){
                    Compound compound = MyCompounds.get(row_model);
                    if(compound!=null){
                        if(!compound.getChemist().equals(userName)){
                            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,"You can only change your own idea");
                            return;
                        }
                        MyCompoundModifyDialog dialog = new MyCompoundModifyDialog(compound, project);
                        dialog.setLocationRelativeTo(MyCompoundIdealPanel.this);
                        dialog.setVisible(true);
                        if(compound.isChanged()){
                            try {
                                calculateProperties(Collections.singletonList(compound));
                                findBioNumbers(Collections.singletonList(compound));
                                return;
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                            }
                        }
                        tableModel.fireTableRowsUpdated(row_model,row_model);
                    }
                }

            }
        });

        commitBtn = new JButton("Save changes.");
        commitBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(compoundsToRemove.size()>0){
                    int option = JOptionPane.showConfirmDialog(MyCompoundIdealPanel.this,"Deleting "+compoundsToRemove.size()+" ideas, are you sure?");
                    if(option!=JOptionPane.YES_OPTION){
                        return;
                    }
                }
                progressMonitor.setProgress(1);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        if(compoundsToRemove.size()>0){
                            Vector v = new Vector();
                            v.add("Deleting ideas ...");
                            v.add(1);
                            publish(v);
                            FrontierDAO.getInstance().deleteMyCompounds(compoundsToRemove);
                        }
                        compoundsToRemove.clear();
                        int idx = 0;
                        for(Compound c:MyCompounds){
                            idx++;
                            if(c.isChanged()){
                                c.save();
                                Vector v = new Vector();
                                v.add("Updating "+ idx +" ideas ...");
                                v.add(100*idx/MyCompounds.size());
                                publish(v);
                            }
                        }
                        return null;
                    }

                    @Override
                    protected void process(List chunks) {
                        Vector v = (Vector) chunks.get(chunks.size() - 1);
                        String note = (String) v.get(0);
                        int progress = (Integer) v.get(1);
                        progressMonitor.setProgress(progress);
                        progressMonitor.setNote(note);
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            tableModel.fireTableDataChanged();
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                        } catch (ExecutionException e2) {
                            e2.printStackTrace();
                            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e2.getMessage());
                        }finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();

            }
        });

        final JCheckBox substructureFilterCB = new JCheckBox("Substructure Filter");
        substructureFilterCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JCheckBox cb = (JCheckBox)e.getSource();
                if(cb.isSelected()){
                    if(substructureSearchDialog==null){
                        substructureSearchDialog = new SubstructureSearchDialog(inSilicoTools.getInstance(),"Substructure",false, false);
                    }
                    substructureSearchDialog.setLocationRelativeTo(MyCompoundIdealPanel.this);
                    substructureSearchDialog.setVisible(true);
                    if(substructureSearchDialog.isCommitted()){
                        String smarts = substructureSearchDialog.getSmarts();
                        if(!smarts.isEmpty()){
                            useSubStructureFilter = true;
                            substructure = smarts;
                            SwingWorker sw = new SwingWorker() {
                                @Override
                                protected Object doInBackground() throws Exception {
                                    substructureMap.clear();
                                    OESubSearch subSearch = new OESubSearch();
                                    subSearch.Init(substructure);
                                    for(Compound c:originalCompounds){
                                        OEGraphMol mol = c.getPropertyMol().getMol();
                                        oechem.OEPrepareSearch(mol,subSearch);
                                        if(subSearch.SingleMatch(mol)){
                                            substructureMap.put(c.getId(),true);
                                        }else{
                                            substructureMap.put(c.getId(),false);
                                        }
                                    }
                                    return null;
                                }

                                @Override
                                protected void done() {
                                    try {
                                        get();
                                        updateTable();
                                    } catch (InterruptedException | ExecutionException e1) {
                                        e1.printStackTrace();
                                        JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                                    }
                                }
                            };
                            sw.execute();
                        }
                    }else{
                        cb.setSelected(false);
                        useSubStructureFilter=false;
                        substructure = null;
                        updateTable();
                    }
                }else{
                    substructure = null;
                    useSubStructureFilter = false;
                    updateTable();
                }

            }
        });

        Vector<Status> My_status = new Stack<Status>();
        My_status.addAll(FrontierDAO.getInstance().getMy_status());
        My_status.insertElementAt(new Status(-1,"All"),0);
        final JComboBox statusFilterCB = new JComboBox(My_status);
        statusFilterCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                statusDesired = (Status) statusFilterCB.getSelectedItem();
                updateTable();
            }
        });


        p.add(addNewBtn);
        p.add(insertBatchBtn);
        p.add(removeBtn);
        p.add(editBtn);
        p.add(commitBtn);
        JPanel filterPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        filterPanel.setMaximumSize(new Dimension(550,50));
        filterPanel.setBorder(new TitledBorder("Filter"));
        filterPanel.add(substructureFilterCB);
        filterPanel.add(new JLabel("Status:"));
        filterPanel.add(statusFilterCB);
        filterPanel.add(new JLabel("Chemist:"));
        filterPanel.add(operationCB);
        filterPanel.add(chemistCB);
        chemistCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                updateTable();
            }
        });
        operationCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                updateTable();
            }
        });
        p.add(filterPanel);

        return p;
    }

    public void updateTable(){
        MyCompounds.clear();
        String chemist = chemistCB.getItemAt(chemistCB.getSelectedIndex());
        boolean showChemist = operationCB.getItemAt(operationCB.getSelectedIndex()).equals("Show");
        for(Compound c:originalCompounds){
            if(showChemist){
                if(!chemist.equals("All")){
                    if(!c.getChemist().equals(chemist)){
                        continue;
                    }
                }
            }else{
                if(chemist.equals("All")){
                    continue;
                }
                if(c.getChemist().equals(chemist)){
                    continue;
                }
            }
            if(statusDesired!=null&&!c.getStatus().equals(statusDesired)&&!statusDesired.getStatus_name().equalsIgnoreCase("All")){
                continue;
            }
            if(useSubStructureFilter&&!substructureMap.get(c.getId())){
                continue;
            }
            MyCompounds.add(c);
        }
        tableModel.fireTableDataChanged();
    }

    public void setProject(Project project) throws SQLException {
        this.project = project;
        reload();
    }

    void reload() throws SQLException {
        final int selectedRow = table.getSelectedRow();
        progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
        SwingWorker sw = new SwingWorker() {
            @Override
            protected Object doInBackground() throws Exception {
                Vector<Compound> MyCompounds = FrontierDAO.getInstance().getMyCompounds(project.getProject_id());
                calculateProperties(MyCompounds);
                return MyCompounds;
            }

            @Override
            protected void done() {
                try{
                    Vector<Compound> compounds = (Vector<Compound>)get();
                    Vector<String> chemists = new Vector<>();
                    for(Compound cmpd:compounds){
                        String chemist = cmpd.getChemist();
                        if(!chemists.contains(chemist)){
                            chemists.add(chemist);
                        }
                    }
                    MyCompounds.clear();
                    MyCompounds.addAll(compounds);
                    originalCompounds.clear();
                    originalCompounds.addAll(compounds);
                    chemistCB.setModel(new DefaultComboBoxModel(chemists));
                    tableModel.fireTableDataChanged();
                    if(selectedRow>=0&&selectedRow<table.getRowCount()) {
                        TableUtils.scrollToVisible(table, selectedRow);
                    }
                    findBioNumbers(compounds);
                } catch (InterruptedException | ExecutionException e) {
                    e.printStackTrace();
                } finally {
                    progressMonitor.close();
                }
            }
        };
        sw.execute();
    }

    class StatusComboboxEditor extends AbstractCellEditor implements TableCellEditor,ActionListener{
        Status status;
        Vector<Status> all_status;
        int row = -1;

        public StatusComboboxEditor(Vector<Status> all_status) {
            this.all_status = all_status;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            JComboBox<Status> statusCB = (JComboBox<Status>) e.getSource();
            this.status = statusCB.getItemAt(statusCB.getSelectedIndex());
            fireEditingStopped();
            tableModel.fireTableRowsUpdated(row,row);
        }

        @Override
        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
            if (value instanceof Status) {
                this.status = (Status) value;
            }
            this.row = row;
            JComboBox<Status> statusCB = new JComboBox<> (all_status);
            statusCB.setSelectedItem(value);
            statusCB.addActionListener(this);

//            if (isSelected) {
//                statusCB.setBackground(table.getSelectionBackground());
//            } else {
//                statusCB.setBackground(table.getSelectionForeground());
//            }

            return statusCB;
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
            return this.status;
        }
    }

    class AssignedChemistComboboxEditor extends AbstractCellEditor implements TableCellEditor,ActionListener{
        Chemist chemist;
        Vector<Chemist> all_chemist;
        int row = -1;

        public AssignedChemistComboboxEditor(Vector<Chemist> all_chemist) {
            this.all_chemist = all_chemist;
        }

        @Override
        public void actionPerformed(ActionEvent e) {
            JComboBox<Chemist> chemistCB = (JComboBox<Chemist>) e.getSource();
            this.chemist = chemistCB.getItemAt(chemistCB.getSelectedIndex());
            fireEditingStopped();
            tableModel.fireTableRowsUpdated(row,row);
        }

        @Override
        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
            if (value instanceof Chemist) {
                this.chemist = (Chemist) value;
            }
            this.row = row;
            if(all_chemist==null){
                all_chemist = new Vector<>();
            }
            JComboBox<Chemist> chemistCB = new JComboBox<> (all_chemist);
            chemistCB.setSelectedItem(value);
            chemistCB.addActionListener(this);
            return chemistCB;
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
            return this.chemist;
        }
    }

    public JMenuBar getJMenuBar(){
        JMenuBar menuBar = new JMenuBar();
        JMenu fileMenu = new JMenu("File");
        JMenuItem changePrjItem = new JMenuItem("Reload from database...");
        changePrjItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    setProject(project);
                } catch (SQLException e1) {
                    e1.printStackTrace();
                    JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                }
            }
        });
        fileMenu.add(changePrjItem);

        JMenuItem exportAllToSDFItem = new JMenuItem("Export all to SDF");
        exportAllToSDFItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveIdeaAsSdf(false);
            }

        });
        fileMenu.add(exportAllToSDFItem);

        JMenuItem exportSelectedToSDFItem = new JMenuItem("Export selected to SDF");
        exportSelectedToSDFItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                saveIdeaAsSdf(true);
            }
        });
        fileMenu.add(exportSelectedToSDFItem);
        menuBar.add(fileMenu);

        JMenu editMenu = new JMenu("Edit");

        JMenuItem selectallItem = new JMenuItem("Select All");
        selectallItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(Compound c:MyCompounds){
                    c.setSelected(true);
                }
            }
        });
        editMenu.add(selectallItem);

        JMenuItem unselectallItem = new JMenuItem("Unselect All");
        unselectallItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(Compound c:MyCompounds){
                    c.setSelected(false);
                }
            }
        });
        editMenu.add(unselectallItem);

        JMenuItem invertselectionItem = new JMenuItem("Invert selection");
        invertselectionItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(Compound c:MyCompounds){
                    if(c.getSelected()){
                        c.setSelected(false);
                    }else{
                        c.setSelected(true);
                    }
                }
            }
        });
        editMenu.add(invertselectionItem);

        JMenuItem selectHighlighedItem = new JMenuItem("Select highlighted");
        selectHighlighedItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(int idx:table.getSelectedRows()){
                    if(idx>=0) {
                        int modelIdx = table.convertRowIndexToModel(idx);
                        if(modelIdx>=0){
                            MyCompounds.get(modelIdx).setSelected(true);
                        }
                    }
                }
            }
        });
        editMenu.add(selectHighlighedItem);

        menuBar.add(editMenu);
        return menuBar;
    }

    void saveIdeaAsSdf(final boolean selectedOnly) {
        if(MyCompounds.size()==0){
            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,"No molecules available.");
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
                        for(Compound c:MyCompounds){
                            if(selectedOnly&&!c.getSelected()){
                                progress ++;
                                continue;
                            }
                            Vector v = new Vector();
                            v.add(String.format("Saving molecule No. %d",progress));
                            v.add(100*progress/MyCompounds.size());
                            publish(v);
                            OEGraphMol mol1 = new OEGraphMol(c.getPropertyMol().getMol());
                            for(int idx=2;idx<tableModel.getColumnCount();idx++){
                                String colName = tableModel.getColumnName(idx);
                                String value = tableModel.getValueAsString(progress, idx);
                                if(value!=null){
                                    oechem.OESetSDData(mol1, colName, value);
                                }else{
                                    oechem.OESetSDData(mol1, colName, "");
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
                            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this, "SDF file saved.");
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(MyCompoundIdealPanel.this,e1.getMessage());
                        }finally {
                            progressMonitor.close();
                        }
                    }
                };
                sw.execute();

            }
        });
    }


    public static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                System.out.println(System.getProperty("os.name"));
                System.out.println(System.getProperty("os.arch"));
                JPopupMenu.setDefaultLightWeightPopupEnabled(false);
                ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);
                if(System.getProperty("os.name").equals("Linux")){
                    try {
                        UIManager.setLookAndFeel(new MetalLookAndFeel());
                    } catch (UnsupportedLookAndFeelException e) {
                        e.printStackTrace();
                    }
                }else {
                    try {
                        PlasticLookAndFeel.setPlasticTheme(new SkyBlue());
                        UIManager.setLookAndFeel(new Plastic3DLookAndFeel());
                    } catch (UnsupportedLookAndFeelException e) {
                        e.printStackTrace();
                    }
                }

                try {
                    OEChemWebLicenseInstaller.loadOELicenseFromWeb();
//                    LicenseManager.setLicenseFile("/Users/jfeng1/.chemaxon/license.cxl");
                    LicenseManager.setLicenseFile("http://10.74.2.128:8080/medchem_design/license.cxl");
                } catch (IOException e) {
                    JOptionPane.showMessageDialog(null,e.getMessage());
                    return;
                } catch (LicenseProcessingException e) {
                    JOptionPane.showMessageDialog(null,e.getMessage());
                    return;
                }
                JFrame f = new JFrame();
                Vector<Compound> cmpds = null;
                try {
                    cmpds = FrontierDAO.getInstance().getMyCompounds(0);
                } catch (SQLException e) {
                    e.printStackTrace();
                    return;
                }
                Compound c = new Compound(0,"0","BIO-0912734\n" +
                        "  Mrv1722707311710292D          \n" +
                        "\n" +
                        " 19 22  0  0  0  0            999 V2000\n" +
                        "    8.4114    8.9057    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    7.0794    7.9342    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    8.4972    8.0784    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    7.6664    9.2388    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    7.8347    7.5978    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    9.1632    9.2318    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    6.4135    7.4502    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    6.9937    8.7547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    6.4993    6.6263    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    9.7159    8.6276    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    9.3039    7.9068    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    5.8333    6.1492    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    5.6582    7.7832    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    7.2512    6.2934    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    4.9923    7.2957    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    5.0815    6.4821    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    5.9191    5.3253    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    7.3403    5.4798    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "    6.6710    4.9923    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                        "  2  8  1  0  0  0  0\n" +
                        "  3  1  1  0  0  0  0\n" +
                        "  4  1  1  0  0  0  0\n" +
                        "  5  3  1  0  0  0  0\n" +
                        "  6  1  2  0  0  0  0\n" +
                        "  7  2  1  0  0  0  0\n" +
                        "  8  4  2  0  0  0  0\n" +
                        "  9  7  1  0  0  0  0\n" +
                        " 10  6  1  0  0  0  0\n" +
                        " 11  3  1  0  0  0  0\n" +
                        " 12  9  2  0  0  0  0\n" +
                        " 13  7  2  0  0  0  0\n" +
                        " 14  9  1  0  0  0  0\n" +
                        " 15 13  1  0  0  0  0\n" +
                        " 16 15  2  0  0  0  0\n" +
                        " 17 12  1  0  0  0  0\n" +
                        " 18 14  2  0  0  0  0\n" +
                        " 19 18  1  0  0  0  0\n" +
                        " 11 10  2  0  0  0  0\n" +
                        "  5  2  2  0  0  0  0\n" +
                        " 16 12  1  0  0  0  0\n" +
                        " 19 17  2  0  0  0  0\n" +
                        "M  END\n","C1=CN2C=C(C=NC2=N1)C1=CC=CC2=C1C=CC=C2",null, new Date(System.currentTimeMillis()),0,0, 0,1);
                cmpds.add(c);
                MyCompoundIdealPanel p = new MyCompoundIdealPanel(new Project(1,"gsk"));
                f.setJMenuBar(p.getJMenuBar());
                f.getContentPane().add(p);
                f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                f.setSize(new Dimension(1280,1024));
                f.setVisible(true);

            }
        });



    }
}
