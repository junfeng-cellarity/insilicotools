package com.insilico.application.insilicotools.gui;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.marvin.beans.MViewPane;
import chemaxon.marvin.io.formats.smiles.SmartsAtomQuerifier;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.PeriodicSystem;
import chemaxon.struc.RxnMolecule;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.data.ReactionRules;
import com.insilico.application.insilicotools.data.ReagentType;
import com.insilico.application.insilicotools.data.SerializableMol;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.gui.table.CompoundPickingPanel;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;
import com.insilico.application.insilicotools.gui.table.MolMatrixTablePane;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.Enumerator;
import com.chemaxon.mapper.AutoMapper;
import com.github.cjwizard.*;
import openeye.oechem.*;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;


/**
 * Created by jfeng1 on 2/2/16.
 */
public class EnumerationWizard extends JPanel {
    EnumerationPageFactory pageFactory;
    WizardContainer wizard;
    WizardSettings settings;
    File currentDirectory;
    Vector<MatrixMolTableModel> unselectedReagentTableModels;
    Vector<MatrixMolTableModel> selectedReagentTableModels;
    Vector<Vector<PropertyMolecule>> unselectedReagentLists;
    Vector<Vector<PropertyMolecule>> selectedReagentLists;

    Vector<PropertyMolecule> products = new Vector<PropertyMolecule>();
    MatrixMolTableModel productTableModel = new MatrixMolTableModel(products);
    Enumerator enumerator;
    RxnPage rxnPage = new RxnPage();
    ReactantPage reactantPage = new ReactantPage();
    ProductPage productPage = new ProductPage();
    DesignProgressMonitor progressMonitor = new DesignProgressMonitor(EnumerationWizard.this,"Progress","Progress",0,100);
    HashMap<String,PropertyMolecule> reagentDb = new HashMap<String, PropertyMolecule>();

    public EnumerationWizard() {
        super(new BorderLayout());
        unselectedReagentLists = new Vector<Vector<PropertyMolecule>>();
        unselectedReagentTableModels = new Vector<MatrixMolTableModel>();
        selectedReagentLists = new Vector<Vector<PropertyMolecule>>();
        selectedReagentTableModels = new Vector<MatrixMolTableModel>();
        enumerator = new Enumerator();
        settings = new StackWizardSettings();
        pageFactory = new EnumerationPageFactory();
        wizard = new WizardContainer(pageFactory){
            @Override
            public void finish() {
                super.finish();
                EnumerationWizard.this.firePropertyChange("Finished",true,false);
            }

            @Override
            public void cancel() {
                super.cancel();
                EnumerationWizard.this.firePropertyChange("Cancelled",true,false);
            }
        };
        wizard.addWizardListener(new WizardListener() {
            @Override
            public void onPageChanged(WizardPage wizardPage, List<WizardPage> list) {
                if(wizardPage instanceof ProductPage){
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker swingWorker = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            String smirks = rxnPage.sketcher.getMol("smarts:ah");
                            enumerator.setSmirks(smirks);
//                            int numReactants = enumerator.getNumReactants();
                            Vector v = new Vector();
                            v.add("Adding reagents");
                            v.add(DesignProgressMonitor.INDETERMINATE);
                            publish(v);
                            reagentDb.clear();
                            for(int i=0;i<selectedReagentLists.size();i++){
                                Vector<PropertyMolecule> reagents = selectedReagentLists.get(i);
                                enumerator.addReagent(reagents,i);
                                for(PropertyMolecule m:reagents){
                                    if(!m.getName().isEmpty()&&!reagentDb.containsKey(m.getName())) {
                                        reagentDb.put(m.getName(), m);
                                    }
                                }
                            }
                            products.clear();
                            v.clear();
                            v.add("Enumerating products");
                            v.add(DesignProgressMonitor.INDETERMINATE);
                            publish(v);
                            Vector<PropertyMolecule> rawProducts = enumerator.enumerate();
//                            v.add("Calculating properties ...");
//                            v.add(DesignProgressMonitor.INDETERMINATE);
//                            publish(v);
//                            ChemFunc.calculateOEProperty(rawProducts);
                            for(PropertyMolecule p: rawProducts){
                                String[] reagentNames = p.getName().split("_");
                                for(int i=0;i<reagentNames.length;i++){
                                    String reagentName = reagentNames[i];
                                    PropertyMolecule r = reagentDb.get(reagentName);
                                    if(r!=null) {
                                        oechem.OESetSDData(p.getMol(), String.format("R_%d_name", i + 1), r.getName());
                                        oechem.OESetSDData(p.getMol(), String.format("R_%d_smiles", i + 1), r.getSmiles());
                                    }
                                }
                                products.add(p);
                            }
                            return products;
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
                            try {
                                get();
                                if(products.size()==0){
                                    throw new Exception("Failed to enumerate compounds, please check the reaction.");
                                }
                                productTableModel.fireTableStructureChanged();
                                progressMonitor.close();
                                JOptionPane.showMessageDialog(EnumerationWizard.this, String.format("%d moelcules enumerated.",products.size()));
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(EnumerationWizard.this,e1.getMessage());
                            }finally {
                                progressMonitor.setNote("Progress");
                                progressMonitor.setProgress(1);
                                progressMonitor.close();
                            }
                        }

                    };
                    swingWorker.execute();


                }else if(wizardPage instanceof ReactantPage){
                    Molecule rxnMol = reactantPage.viewPane.getM(0);
                    if(rxnMol.isReaction()){
                        RxnMolecule mol = RxnMolecule.getReaction(rxnMol);
                        mol.setAtomSetSeqs(0);
                        Molecule reactant = mol.getReactant(reactantPage.tabbedPane.getSelectedIndex());
                        reactant.setAtomSetSeqs(1);
                        reactantPage.viewPane.setM(0,mol);
                    }

                }
            }

            @Override
            public void onFinished(List<WizardPage> list, WizardSettings wizardSettings) {
            }

            @Override
            public void onCanceled(List<WizardPage> list, WizardSettings wizardSettings) {
            }
        });
        add(wizard,BorderLayout.CENTER);
    }

    public void loadSession(String smirks, ArrayList<ArrayList<ArrayList<SerializableMol>>> reagents, ArrayList<SerializableMol> productList) throws MolFormatException {
        enumerator.setSmirks(smirks);
        Molecule rxnMol = MolImporter.importMol(smirks);
        rxnPage.sketcher.setMol(rxnMol);
        rxnPage.sketcher.firePropertyChange("mol",true,false);
//        reactantPage.setReactant(enumerator.getNumReactants(), rxnMol);
        products.clear();
        boolean reagentsSelected = true;
        if(smirks!=null&&smirks.length()>0&&reagents.size()>0){
            for(int i=0;i<reagents.size();i++){
                ArrayList<ArrayList<SerializableMol>> reactants = reagents.get(i);
                Vector<PropertyMolecule> selectedList = new Vector<PropertyMolecule>();
                Vector<PropertyMolecule> unselectedList = new Vector<PropertyMolecule>();
                for(int j = 0; j< reactants.size(); j++){
                    ArrayList<SerializableMol> unSelectedReactants = reactants.get(0);
                    for(SerializableMol m:unSelectedReactants){
                        unselectedList.add(new PropertyMolecule(m.getOEMol()));
                    }
                    ArrayList<SerializableMol> selectedReactants = reactants.get(1);
                    for(SerializableMol m:selectedReactants){
                        selectedList.add(new PropertyMolecule(m.getOEMol()));
                    }
                    if(selectedReactants.size()==0&& reagentsSelected){
                        reagentsSelected = false;
                    }
                }
                selectedReagentLists.get(i).addAll(selectedList);
                selectedReagentTableModels.get(i).fireTableStructureChanged();
                unselectedReagentLists.get(i).addAll(unselectedList);
                unselectedReagentTableModels.get(i).fireTableDataChanged();
            }
            reactantPage.setReagentsLoaded(reagentsSelected);
            for(SerializableMol m:productList){
                products.add(new PropertyMolecule(m.getOEMol()));
            }
            productTableModel.fireTableDataChanged();
        }
    }

    public Vector<PropertyMolecule> getProducts() {
        return products;
    }

    private void loadSdfToTable(final Vector<PropertyMolecule> molecules, final MatrixMolTableModel tableModel){
        JFileChooser fc = new JFileChooser(currentDirectory);
        fc.setFileFilter(new FileNameExtensionFilter("SDF file.", "sdf"));
        int option = fc.showOpenDialog(EnumerationWizard.this);
        if (option == JFileChooser.APPROVE_OPTION) {
            currentDirectory = fc.getCurrentDirectory();
            File f = fc.getSelectedFile();
            if (f.canRead()) {
                molecules.clear();
                final String fname = f.getAbsolutePath();
//                    OEMolDatabase db = new OEMolDatabase(fname);
//                    final int numMols = db.GetMaxMolIdx();
                progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        oemolistream ifs = new oemolistream();
                        ifs.SetFormat(OEFormat.SDF);
                        ifs.open(fname);
                        Vector<String> properties = new Vector<String>();
                        OEGraphMol mol = new OEGraphMol();
                        while (oechem.OEReadMolecule(ifs, mol)) {
                            OESDDataIter iter = oechem.OEGetSDDataPairs(mol);
                            while(iter.hasNext()){
                                OESDDataPair pair = iter.next();
                                if(!properties.contains(pair.GetTag())) {
                                    properties.add(pair.GetTag());
                                }
                            }
                            molecules.add(new PropertyMolecule(mol));
                        }
//                        ChemFunc.calculateOEProperty(molecules);
//                        ChemFunc.generateCLogP(molecules);
//
                        return properties;
                    }

                    @Override
                    protected void done() {
                        try {
                            Vector<String> properties = (Vector<String>)get();
                            properties.insertElementAt("No Name",0);
                            String[] tags = properties.toArray(new String[properties.size()]);
                            progressMonitor.close();
                            String selectedTag = (String) JOptionPane.showInputDialog(EnumerationWizard.this, "Input", "Mol Name", JOptionPane.QUESTION_MESSAGE, null, tags, tags[0]);
                            for (PropertyMolecule m : molecules) {
                                if(selectedTag.equals("No Name")){
                                    m.setName("");
                                }else {
                                    String title = oechem.OEGetSDData(m.getMol(), selectedTag);
                                    m.setName(title);
                                }
                            }
                            tableModel.fireTableStructureChanged();
                        } catch (InterruptedException e1) {
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(EnumerationWizard.this, e1.getMessage());
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            progressMonitor.close();
                            JOptionPane.showMessageDialog(EnumerationWizard.this, e1.getMessage());
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

    private void saveSession() {
        FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("Library Session file", "enum"), new FileFunctor() {
            @Override
            public void execute(final File file) {
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        String smirks = enumerator.getSmirks();
                        ArrayList<ArrayList<ArrayList<SerializableMol>>> reagents = new ArrayList<ArrayList<ArrayList<SerializableMol>>>();
                        for(int i=0;i<enumerator.getNumReactants();i++){
                            ArrayList<ArrayList<SerializableMol>> reactants = new ArrayList<ArrayList<SerializableMol>>();
                            Vector<PropertyMolecule> p1 = unselectedReagentLists.get(i);
                            Vector<PropertyMolecule> p2 = selectedReagentLists.get(i);
                            ArrayList<SerializableMol> p1s = new ArrayList<SerializableMol>();
                            for(PropertyMolecule p:p1){
                                p1s.add(p.getSerializableMol2D());
                            }
                            ArrayList<SerializableMol> p2s = new ArrayList<SerializableMol>();
                            for(PropertyMolecule p:p2){
                                p2s.add(p.getSerializableMol2D());
                            }
                            reactants.add(p1s);
                            reactants.add(p2s);
                            reagents.add(reactants);
                        }
                        String outputfile = file.getAbsolutePath();
                        if(!outputfile.endsWith(".enum")){
                            outputfile = outputfile.split("\\.")[0]+".enum";
                        }
                        ArrayList<SerializableMol> productList = new ArrayList<SerializableMol>();
                        int progress = 0;
                        for (PropertyMolecule mol : products) {
                            progress++;
                            Vector v = new Vector();
                            v.add(String.format("Saving molecule No. %d", progress));
                            v.add(100 * progress / products.size());
                            publish(v);
                            productList.add(mol.getSerializableMol2D());
                        }
                        FileOutputStream fos = new FileOutputStream(outputfile);
                        ObjectOutputStream oos = new ObjectOutputStream(fos);
                        oos.writeObject(smirks);
                        oos.writeObject(reagents);
                        oos.writeObject(productList);
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
                            JOptionPane.showMessageDialog(EnumerationWizard.this, String.format("Session file saved as %s.",output));
                        } catch (Exception e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(EnumerationWizard.this, e1.getMessage());
                        } finally {
                        }
                    }
                };
                sw.execute();

            }
        });
    }


    private class ProductPage extends WizardPage{

        public ProductPage() {
            super("Products", "Enumerated Products");
            setLayout(new BorderLayout());
            MolMatrixTablePane matrixTablePane = new MolMatrixTablePane(productTableModel);
            JPanel btnPanel = matrixTablePane.getBtnPanel();
            JButton saveAsSessionBtn = new JButton("Save As Session");
            saveAsSessionBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    saveSession();
                }
            });
            btnPanel.add(saveAsSessionBtn);
            add(matrixTablePane,BorderLayout.CENTER);
        }

        @Override
        public void rendering(List<WizardPage> path, WizardSettings settings) {
            super.rendering(path, settings);
            setNextEnabled(false);
            setPrevEnabled(true);
            setFinishEnabled(true);
        }
    }

    private class ReactantPage extends WizardPage{
        int numReactants;
        MViewPane viewPane;
        boolean isReagentsLoaded = false;
        JTabbedPane tabbedPane;

        public void setReagentsLoaded(boolean reagentsLoaded) {
            isReagentsLoaded = reagentsLoaded;
        }

        public ReactantPage() {
            super("Reactants", "Load Reactants ...");
            viewPane = MarvinFactory.createViewPane();
            viewPane.setEditable(MViewPane.VIEW_ONLY);
            setLayout(new BorderLayout());
        }

        @Override
        public void rendering(List<WizardPage> path, WizardSettings settings) {
            super.rendering(path, settings);
            setNextEnabled(isReagentsLoaded);
            setPrevEnabled(true);
            setFinishEnabled(false);
        }

        public void setReactant(int numReactants, Molecule rxnMol) {
            removeAll();
            this.numReactants = numReactants;
            selectedReagentLists.clear();
            selectedReagentTableModels.clear();
            unselectedReagentLists.clear();
            unselectedReagentTableModels.clear();
            if(numReactants==0){
                add(new JLabel("No reaction available, please redraw the reaction."));
            }else{
                setLayout(new BorderLayout());
                add(viewPane,BorderLayout.NORTH);
                viewPane.setPreferredSize(new Dimension(1000,100));
                viewPane.setM(0,rxnMol);
                tabbedPane = new JTabbedPane();
                tabbedPane.addChangeListener(new ChangeListener() {
                    @Override
                    public void stateChanged(ChangeEvent e) {
                        Molecule rxnMol = viewPane.getM(0);
                        if(rxnMol.isReaction()){
                            RxnMolecule mol = RxnMolecule.getReaction(rxnMol);
                            mol.setAtomSetSeqs(0);
                            Molecule reactant = mol.getReactant(tabbedPane.getSelectedIndex());
                            reactant.setAtomSetSeqs(1);
                            viewPane.setM(0,mol);
                        }
                    }
                });
                for(int i=0;i<numReactants;i++){
                    JPanel pi = new JPanel(new BorderLayout());
                    final Vector<PropertyMolecule> unselected_reagents = new Vector<PropertyMolecule>();
                    final Vector<PropertyMolecule> selected_reagents = new Vector<PropertyMolecule>();
                    selectedReagentLists.add(selected_reagents);
                    unselectedReagentLists.add(unselected_reagents);
                    CompoundPickingPanel reagentPanel = new CompoundPickingPanel(unselected_reagents,selected_reagents);
                    final MatrixMolTableModel leftTableModel = reagentPanel.getLeftTableModel();
                    final MatrixMolTableModel rightTableModel = reagentPanel.getRightTableModel();
                    final int tabIdx = i;
                    leftTableModel.addTableModelListener(new TableModelListener() {
                        @Override
                        public void tableChanged(TableModelEvent e) {
                            int rightSize = rightTableModel.getPropertyMolecules().size();
                            int leftSize = leftTableModel.getPropertyMolecules().size();
                            tabbedPane.setTitleAt(tabIdx,String.format("Reagents %d (%d out of %d selected)",tabIdx, rightSize, rightSize+leftSize));
                        }
                    });

                    rightTableModel.addTableModelListener(new TableModelListener() {
                        @Override
                        public void tableChanged(TableModelEvent e) {
                            Vector<PropertyMolecule> propertyMolecules = rightTableModel.getPropertyMolecules();
                            for(PropertyMolecule m:propertyMolecules){
                                if(!selectedReagentLists.get(tabIdx).contains(m)){
                                    selectedReagentLists.get(tabIdx).add(m);
                                }
                            }
                            boolean allReagentsLoaded = true;
                            for(Vector<PropertyMolecule> reagents1:selectedReagentLists){
                                if(reagents1.size()==0){
                                    allReagentsLoaded = false;
                                    break;
                                }
                            }
                            if(allReagentsLoaded){
                                isReagentsLoaded = true;
                                setNextEnabled(true);
                            }else{
                                isReagentsLoaded = false;
                                setNextEnabled(false);
                            }
                        }
                    });

                    unselectedReagentTableModels.add(leftTableModel);
                    selectedReagentTableModels.add(rightTableModel);
                    pi.add(reagentPanel,BorderLayout.CENTER);
                    tabbedPane.addTab(String.format("Reagents %d",i),pi);
                }

                JPanel btnPanel = new JPanel();
                JButton searchReagentBtn = new JButton("Search Frontier");
                searchReagentBtn.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        Container topLevelAncestor = EnumerationWizard.this.getTopLevelAncestor();
                        final SubstructureSearchDialog dialog;
                        if(topLevelAncestor instanceof JDialog){
                            dialog = new SubstructureSearchDialog((JDialog)topLevelAncestor,"Search Frontier Reagent Catalogue");
                        }else{
                            dialog = new SubstructureSearchDialog((JFrame)topLevelAncestor,"Search Frontier Reagent Catalogue");

                        }
                        dialog.setLocationRelativeTo(EnumerationWizard.this);
                        dialog.setVisible(true);
                        final String smarts = dialog.getSmarts();
                        final int hits = dialog.getMaxNumHits();
                        final float clogp = dialog.getCLogPLimit();
                        final float mw = dialog.getMWLimit();
                        final float psa = dialog.getPSALimit();

                        if(dialog.isCommitted()&&!smarts.isEmpty()){
                            progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                            SwingWorker sw = new SwingWorker() {
                                @Override
                                protected Object doInBackground() throws Exception {
                                    Vector<String> molFromSmarts = FrontierDAO.getInstance().getMolFromSmartsFromFrontier(smarts,hits,clogp,mw,psa);
                                    Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
                                    for(String molString:molFromSmarts){
                                        oemolistream ifs = new oemolistream();
                                        ifs.SetFormat(OEFormat.SDF);
                                        ifs.openstring(molString);
                                        OEGraphMol oemol = new OEGraphMol();
                                        oechem.OEReadMolecule(ifs, oemol);
                                        propertyMolecules.add(new PropertyMolecule(oemol));
                                    }
                                    ChemFunc.calculateOEProperty(propertyMolecules);
                                    return propertyMolecules;
                                }

                                @Override
                                protected void done() {
                                    try {
                                        Vector<PropertyMolecule> mols = (Vector<PropertyMolecule>)get();
                                        int id = tabbedPane.getSelectedIndex();
                                        Vector<PropertyMolecule> propertyMolecules = unselectedReagentLists.get(id);
                                        MatrixMolTableModel matrixMolTableModel = unselectedReagentTableModels.get(id);
                                        propertyMolecules.addAll(mols);
                                        matrixMolTableModel.fireTableDataChanged();
                                    } catch (Exception e1) {
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
                btnPanel.add(searchReagentBtn);
                JButton loadReagentsBtn = new JButton("Load SDF");
                loadReagentsBtn.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        int id = tabbedPane.getSelectedIndex();
                        Vector<PropertyMolecule> propertyMolecules = unselectedReagentLists.get(id);
                        MatrixMolTableModel matrixMolTableModel = unselectedReagentTableModels.get(id);
                        loadSdfToTable(propertyMolecules,matrixMolTableModel);
                        unselectedReagentLists.get(id).clear();
                        unselectedReagentTableModels.get(id).fireTableDataChanged();
                    }
                });
                btnPanel.add(loadReagentsBtn);

                JButton sketchBtn = new JButton("Sketch");
                sketchBtn.addActionListener(new ActionListener() {
                    @Override
                    public void actionPerformed(ActionEvent e) {
                        int id = tabbedPane.getSelectedIndex();
                        Vector<PropertyMolecule> propertyMolecules = unselectedReagentLists.get(id);
                        MatrixMolTableModel matrixMolTableModel = unselectedReagentTableModels.get(id);
                        CompoundInputDialog dialog = new CompoundInputDialog();
                        dialog.setLocationRelativeTo(EnumerationWizard.this);
                        dialog.setVisible(true);
                        if(dialog.isCommitted){
                            try {
                                String qmolString = MolExporter.exportToFormat(dialog.getMolecule(), "sdf");
                                OEGraphMol qmol = ChemFunc.getMolFromMolString(qmolString,OEFormat.SDF);
                                propertyMolecules.add(new PropertyMolecule(qmol));
                                matrixMolTableModel.fireTableDataChanged();
                            } catch (IOException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(EnumerationWizard.this,e1.getMessage());
                            }
                        }

                    }
                });
                btnPanel.add(sketchBtn);
                add(tabbedPane,BorderLayout.CENTER);
                add(btnPanel,BorderLayout.SOUTH);
            }
        }
    }

    private class RxnPage extends WizardPage{
        MSketchPane sketcher;
        public RxnPage() {
            super("Reaction", "Sketch the reaction for this enumeration");
            sketcher = MarvinFactory.createCompoundSketcher();
            sketcher.addPropertyChangeListener(new PropertyChangeListener() {
                @Override
                public void propertyChange(PropertyChangeEvent evt) {
                    Molecule mol = sketcher.getMol();
                    if(evt.getPropertyName().equals("mol")){
                        if(!mol.isEmpty()&& mol.isReaction()){
                            setNextEnabled(true);
                            if(reactantPage!=null){
                                RxnMolecule rxnMol = RxnMolecule.getReaction(rxnPage.sketcher.getMol());
                                Molecule[] reactants = rxnMol.getReactants();
                                if(reactants.length!=reactantPage.numReactants) {
                                    reactantPage.setReactant(reactants.length, rxnMol);
                                }
                                reactantPage.viewPane.setM(0, mol);
                            }
                        }else{
                            setNextEnabled(false);
                        }
                    }else if(evt.getPropertyName().equals("selectionChanged")){
                        if(!mol.isEmpty()){
                            MolAtom selectedAtm = null;
                            int numSelectedAtm = 0;
                            for(MolAtom atm:mol.getAtomArray()){
                                if(atm.isSelected()){
                                    numSelectedAtm += 1;
                                    selectedAtm = atm;
                                }
                            }
                            if(numSelectedAtm==1){
                                ReactionRules[] options = null;
                                if(selectedAtm.getAtno()== PeriodicSystem.N){
                                    options = ReactionRules.getNitrogenRules();
                                }else if(selectedAtm.getAtno()==PeriodicSystem.O){
                                    options = ReactionRules.getOxygenRules();
                                }else if(selectedAtm.getAtno()==PeriodicSystem.C){
                                    options = ReactionRules.getCarbonRules();
                                }
                                if(options!=null){
                                    ReactionRules rule = (ReactionRules) JOptionPane.showInputDialog(EnumerationWizard.this, "Selected atom is ", "Atom Type", JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
                                    if(rule!=null){
                                        if(!rule.isAny()) {
                                            int mapIdx = selectedAtm.getAtomMap();
                                            SmartsAtomQuerifier.setSMARTS(selectedAtm, rule.getSmarts());
                                            selectedAtm.setAtomMap(mapIdx);
                                            rxnPage.sketcher.setMol(mol);
                                            rxnPage.sketcher.firePropertyChange("mol",true,false);
                                        }else{
                                            int mapIdx = selectedAtm.getAtomMap();
                                            selectedAtm.clearQProps();
                                            selectedAtm.setAtomMap(mapIdx);
                                            rxnPage.sketcher.setMol(mol);
                                            rxnPage.sketcher.firePropertyChange("mol",true,false);
                                        }
                                    }
                                }

                            }
                        }
                    }

                }
            });
            setLayout(new BorderLayout());
            add(sketcher,BorderLayout.CENTER);
            add(buildTopPanel(),BorderLayout.NORTH);
            add(buildBottomPanel(),BorderLayout.SOUTH);
        }

        JPanel buildTopPanel(){
            JPanel p = new JPanel(new FlowLayout(FlowLayout.CENTER,20,-10));
            p.setBorder(new EtchedBorder());
            p.setPreferredSize(new Dimension(500,50));
            JLabel helpTxtLbl = new JLabel("<html><p><h3>Sketch a reaction and use the arrow to map the unchanging atoms.</h3></p></html>");
            helpTxtLbl.setVerticalAlignment(JLabel.CENTER);
            helpTxtLbl.setHorizontalAlignment(JLabel.CENTER);
            p.add(helpTxtLbl);
            return p;
        }

        JPanel buildBottomPanel(){
            JPanel p = new JPanel();
            p.setBorder(new EtchedBorder());
            p.setPreferredSize(new Dimension(500,50));
            JButton mapBtn = new JButton("AutoMap");
            mapBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Molecule mol = sketcher.getMol();
                    RxnMolecule reaction = RxnMolecule.getReaction(mol);
                    if(!mol.isEmpty()){
                        AutoMapper mapper = new AutoMapper();
                        mapper.setMappingStyle(AutoMapper.MappingStyle.MATCHING);
                        mapper.map(reaction);
                        sketcher.setMol(reaction);
                    }
                }
            });
            p.add(mapBtn);

            JButton unmapBtn = new JButton("Remove Atom Mapping");
            unmapBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Molecule mol = sketcher.getMol();
                    if(!mol.isEmpty()){
                        AutoMapper.unmap(mol);
                        sketcher.setMol(mol);
                    }
                }
            });
            p.add(unmapBtn);

            final LinkedHashMap<String,String> reagentsMap = new LinkedHashMap<String, String>();
            reagentsMap.put("Primary Alkyl Amine", ReagentType.PRIMARY_ALKYL_AMINE);
            reagentsMap.put("Primary Amine", ReagentType.PRIMARY_AMINE);
            reagentsMap.put("Secondary Alkyl Amine",ReagentType.SECONDARY_ALKYL_AMINE);
            reagentsMap.put("Secondary Amine",ReagentType.SECONDARY_AMINE);
            reagentsMap.put("Primary/Secondary Alkyl Amine",ReagentType.PRIMARY_SECONDARY_ALKYL_AMINE);
            reagentsMap.put("Primary/Secondary Amine",ReagentType.PRIMARY_SECONDARY_AMINE);
            reagentsMap.put("Aldehyde",ReagentType.ALDEHYDE);
            reagentsMap.put("Carboxylic Acid",ReagentType.CARBOXYLIC_ACID);
            reagentsMap.put("Sulfonyl Chloride",ReagentType.SULFONYL_CHLORIDE);
            reagentsMap.put("Aryl Halide",ReagentType.ARYL_HALIDE);
            reagentsMap.put("Alkyl Halide",ReagentType.ALKYL_HALIDE);
            reagentsMap.put("Halide", ReagentType.HALIDE);
            reagentsMap.put("Chloroformate",ReagentType.CHLOROFORMATE);
            reagentsMap.put("Boronic Acids and Esters",ReagentType.BORONIC_ACIDS_AND_ESTERS);


            final JComboBox reagentListCB = new JComboBox(new Vector<String>(reagentsMap.keySet()));
            reagentListCB.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                    String reagentName = (String)reagentListCB.getSelectedItem();
                    clipboard.setContents(new StringSelection(reagentsMap.get(reagentName)),null);
                    sketcher.doPaste();
                }
            });

            JButton addReagentBtn = new JButton("Add Reagent:");
            addReagentBtn.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                    String reagentName = (String)reagentListCB.getSelectedItem();
                    clipboard.setContents(new StringSelection(reagentsMap.get(reagentName)),null);
                    sketcher.doPaste();
                }
            });
            JPanel p1 = new JPanel();
            p1.setBorder(new LineBorder(Color.BLACK,1));
            p1.add(addReagentBtn);
            p1.add(reagentListCB);
            p.add(p1);
            return p;
        }

        @Override
        public String getTitle() {
            return super.getTitle();
        }

        @Override
        public String getDescription() {
            return super.getDescription();
        }

        @Override
        public void updateSettings(WizardSettings settings) {
            settings.put("reaction", sketcher.getMol());
        }

        @Override
        public boolean onNext(WizardSettings settings) {
            updateSettings(settings);
            return true;
        }

        @Override
        public void rendering(List<WizardPage> path, WizardSettings settings) {
            if(!sketcher.getMol().isEmpty()&&sketcher.getMol().isReaction()){
                setNextEnabled(true);
            }else{
                setNextEnabled(false);
            }
            setPrevEnabled(false);
            setFinishEnabled(false);
        }

    }

    private class EnumerationPageFactory implements PageFactory{
        public EnumerationPageFactory() {

        }

        @Override
        public WizardPage createPage(List<WizardPage> list, WizardSettings wizardSettings) {
            switch (list.size()){
                case 0:
                    return rxnPage;
                case 1:
//                    RxnMolecule rxnMol = RxnMolecule.getReaction(rxnPage.sketcher.getMol());
//                    reactantPage.setReactant(rxnMol.getReactantCount(), rxnMol);
                    return reactantPage;
                case 2:
                    return productPage;
                default:
                    return null;
            }
        }

        @Override
        public boolean isTransient(List<WizardPage> list, WizardSettings wizardSettings) {
            return false;
        }

    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        EnumerationWizard wizard = new EnumerationWizard();
        f.getContentPane().add(wizard);
        f.setSize(new Dimension(1024,768));
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
    }
}
