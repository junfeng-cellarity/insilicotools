package com.insilico.application.insilicotools.gui;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.beans.MViewPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.Core;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import com.insilico.application.insilicotools.gui.widget.CoreStatus;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.google.common.base.Strings;
import com.jidesoft.swing.JideSplitPane;
import com.jidesoft.swing.JideTitledBorder;
import openeye.oechem.OEAtomBase;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OEIsRGroup;
import org.jdesktop.swingx.JXImageView;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by jfeng1 on 2/9/17.
 */
public class MyCoreInputPanel extends JPanel {
    CompoundInputDialog dialog = new CompoundInputDialog();
    CompoundInputDialog reactDialog = new CompoundInputDialog();
    Molecule coreMol;
    Molecule reactionMol;
    String coreSmiles;
    byte[] cdx;
    String description;
    String name;
    int number_of_diversities;
    String coreMolStr;
    JButton saveBtn;
    Core core;
    int status_id;
    JComboBox statusCB;

    public MyCoreInputPanel(){
        this(null);
    }

    public Core getCore() {
        return core;
    }

    public MyCoreInputPanel(Core core1) {
        super(new BorderLayout());
        this.core = core1;
        JPanel namePanel = new JPanel(new FlowLayout(FlowLayout.LEADING));
        namePanel.setBorder(new JideTitledBorder("Name"));
        final JTextField textfield = new JTextField(40);
        if(core!=null){
            name = core.getName();
            textfield.setText(core.getName());
        }
        textfield.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                name = textfield.getText().trim();
                checkAndUpdateStatus();
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                name = textfield.getText().trim();
                checkAndUpdateStatus();
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                name = textfield.getText().trim();
                checkAndUpdateStatus();
            }
        });
        namePanel.add(textfield);

        JideSplitPane p = new JideSplitPane(JideSplitPane.VERTICAL_SPLIT);
        JPanel buttonPanel = new JPanel();
        saveBtn = new JButton("Save Library!");
        saveBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(core==null) {
                    final String chemist = System.getProperty("user.name");
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            int core_id = FrontierDAO.getInstance().insertMyCoreWrapper(name, coreMolStr, coreSmiles, number_of_diversities, description, cdx, chemist, "Cellarity", status_id);
                            return FrontierDAO.getInstance().getCore(core_id);
                        }

                        @Override
                        protected void done() {
                            try {
                                core = (Core)get();
                                MyCoreInputPanel.this.firePropertyChange("LibraryUpdated", true, false);
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                            } catch (ExecutionException e1) {
                                String errorPattern = ".*Detail: Key \\(core_smiles\\)=\\((.*)\\).*$";
                                Pattern p = Pattern.compile(errorPattern);
                                Matcher m = p.matcher(e1.getMessage().replace("\n",""));
                                if(m.matches()){
                                    try {
                                        String smiles = m.group(1);
                                        JDialog dialog = new JDialog();
                                        dialog.setTitle("Library Match Found!");
                                        MViewPane viewPane = MarvinFactory.createViewPane();
                                        viewPane.setM(0,MolImporter.importMol(smiles));
                                        dialog.setModal(true);
                                        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                                        dialog.getContentPane().add(viewPane);
                                        dialog.setSize(new Dimension(800,600));
                                        dialog.setLocationRelativeTo(MyCoreInputPanel.this);
                                        dialog.setVisible(true);
                                    } catch (MolFormatException e2) {
                                        JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                                    }
                                }else{
                                    System.out.println("No Match!");
                                    System.out.println(e1.getMessage());
                                    JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                                }
                            }
                        }
                    };
                    sw.execute();
                }else{
                    final int core_id = core.getId();
                    core.setName(name);
                    core.setCdx(cdx);
                    core.setCore_mol(coreMolStr);
                    core.setCore_smiles(coreSmiles);
                    core.setDescription(description);
                    core.setStatus_id(status_id);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            FrontierDAO.getInstance().updateMyCoreWrapper(core_id,name,coreMolStr,coreSmiles,number_of_diversities,description,cdx,core.getChemist(),core.getSource(),status_id);
                            return null;
                        }
                        @Override
                        protected void done() {
                            try {
                                get();
                                MyCoreInputPanel.this.firePropertyChange("LibraryUpdated", true, false);
                            } catch (InterruptedException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                            } catch (ExecutionException e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                            }
                        }
                    };
                    sw.execute();
                }
            }
        });
        buttonPanel.add(saveBtn);
        saveBtn.setEnabled(false);
        JPanel corePanel = buildCoreImagePanel();
        JPanel reactionPanel = buildReactionImagePanel();
        JPanel descriptionPanel = buildDescriptionPanel();
        JPanel statusPanel = buildStatusPanel();
        p.add(namePanel);
        p.add(corePanel);
        p.add(reactionPanel);
        p.add(descriptionPanel);
        p.add(statusPanel);

        p.setProportionalLayout(true);
        p.setProportions(new double[]{0.1,0.4,0.2,0.2});
        add(p,BorderLayout.CENTER);
        add(buttonPanel,BorderLayout.SOUTH);
    }

    private JPanel buildStatusPanel(){
        JPanel p = new JPanel();
        HashMap<Integer,CoreStatus> statusDict = new HashMap<Integer, CoreStatus>();
        Vector<CoreStatus> allCoreStatus = new Vector<CoreStatus>();
        try {
            allCoreStatus = FrontierDAO.getInstance().getAllCoreStatus();
            for(CoreStatus status:allCoreStatus){
                statusDict.put(status.getStatus_id(),status);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
        statusCB = new JComboBox(allCoreStatus);
        statusCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                status_id = ((CoreStatus)statusCB.getSelectedItem()).getStatus_id();
                checkAndUpdateStatus();
            }
        });
        if(core!=null){
            statusCB.setSelectedItem(statusDict.get(core.getStatus_id()));
            status_id = core.getStatus_id();
        }else{
            status_id = ((CoreStatus)statusCB.getSelectedItem()).getStatus_id();
        }
        p.add(statusCB);
        return p;
    }

    private JPanel buildCoreImagePanel(){
        JPanel p = new JPanel(new BorderLayout());
        final JXImageView imageView = new JXImageView();
        if(core!=null){
            OEGraphMol oemol = core.getCore_mol().getMol();
            coreMol = OEChemFunc.getInstance().convertOEChemMol(oemol);
            coreSmiles = core.getCore_smiles();
            coreMolStr = ChemFunc.getMolString(oemol);
            number_of_diversities = 0;
            for(OEAtomBase atm:oemol.GetAtoms()){
                if(new OEIsRGroup(0).constCall(atm)){
                    number_of_diversities += 1;
                }
            }
            try {
                byte[] png = MolExporter.exportToBinFormat(coreMol, "png");
                imageView.setImage(ImageIO.read(new ByteArrayInputStream(png)));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        JPanel btnPanel = new JPanel();
        JButton addBtn = new JButton("Add/Change Library Core");
        addBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(coreMol!=null){
                    OEGraphMol oemol = OEChemFunc.getInstance().convertChemAxonMol(coreMol);
                    dialog.setMolecule(new PropertyMolecule(oemol));
                }
                dialog.setLocationRelativeTo(MyCoreInputPanel.this);
                dialog.setVisible(true);
                if(dialog.isCommitted()){
                    Molecule mol = dialog.getMolecule();
                    try {
                        byte[] png = MolExporter.exportToBinFormat(mol, "png");
                        imageView.setImage(ImageIO.read(new ByteArrayInputStream(png)));
                        coreMol = mol;
                        OEGraphMol oemol = OEChemFunc.getInstance().convertChemAxonMol(coreMol);
                        number_of_diversities = 0;
                        for(OEAtomBase atm:oemol.GetAtoms()){
                            if(new OEIsRGroup(0).constCall(atm)){
                                number_of_diversities += 1;
                            }
                        }
                        coreSmiles = MolExporter.exportToFormat(coreMol,"smiles:u,a");
                        coreMolStr = ChemFunc.getMolString(oemol);
                        checkAndUpdateStatus();
                    } catch (IOException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                    }
                }
            }
        });
        btnPanel.add(addBtn);
        p.add(imageView,BorderLayout.CENTER);
        p.add(btnPanel,BorderLayout.SOUTH);
        p.setBorder(new TitledBorder("Core"));
        return p;
    }

    private JPanel buildReactionImagePanel(){
        JPanel p = new JPanel(new BorderLayout());
        final JXImageView imageView = new JXImageView();
        if(core!=null){
            cdx = core.getCdx();
            try {
                Molecule mol = MolImporter.importMol(cdx);
                byte[] png = MolExporter.exportToBinFormat(mol, "png");
                imageView.setImage(ImageIO.read(new ByteArrayInputStream(png)));
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        JPanel btnPanel = new JPanel();
        JButton addBtn = new JButton("Add/Change Reaction Scheme");
        addBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(cdx!=null){
                    try {
                        Molecule mol = MolImporter.importMol(cdx);
                        reactDialog.setChemAxonMol(mol);
                    } catch (MolFormatException e1) {
                        e1.printStackTrace();
                    }
                }
                reactDialog.setLocationRelativeTo(MyCoreInputPanel.this);
                reactDialog.setVisible(true);
                if(reactDialog.isCommitted()){
                    Molecule mol = reactDialog.getMolecule();
                    try {
                        byte[] png = MolExporter.exportToBinFormat(mol, "png");
                        imageView.setImage(ImageIO.read(new ByteArrayInputStream(png)));
                        reactionMol = mol;
                        cdx = MolExporter.exportToBinFormat(mol, "cdx");
                        checkAndUpdateStatus();
                    } catch (IOException e1) {
                        e1.printStackTrace();
                        JOptionPane.showMessageDialog(MyCoreInputPanel.this,e1.getMessage());
                    }
                }
            }
        });
        btnPanel.add(addBtn);
        p.add(imageView,BorderLayout.CENTER);
        p.add(btnPanel,BorderLayout.SOUTH);
        p.setBorder(new TitledBorder("Reaction"));
        return p;
    }

    private JPanel buildDescriptionPanel(){
        JPanel p = new JPanel(new BorderLayout());
        final JEditorPane panel = new JEditorPane();
        if(core!=null){
            panel.setText(core.getDescription());
        }
        panel.getDocument().addDocumentListener(new DocumentListener() {
            @Override
            public void insertUpdate(DocumentEvent e) {
                description = panel.getText();
                checkAndUpdateStatus();
            }

            @Override
            public void removeUpdate(DocumentEvent e) {
                description = panel.getText();
                checkAndUpdateStatus();
            }

            @Override
            public void changedUpdate(DocumentEvent e) {
                description = panel.getText();
                checkAndUpdateStatus();
            }
        });
        p.add(panel,BorderLayout.CENTER);
        p.setBorder(new TitledBorder("Description"));
        return p;
    }

    private void checkAndUpdateStatus(){
        if(saveBtn!=null){
            if(Strings.isNullOrEmpty(name)){
                saveBtn.setEnabled(false);
                return;
            }
            if(Strings.isNullOrEmpty(coreMolStr)){
                saveBtn.setEnabled(false);
                return;
            }
            if(cdx==null||cdx.length==0){
                saveBtn.setEnabled(false);
                return;
            }
            saveBtn.setEnabled(true);
        }
    }

    public static void main(String[] args) {
        String username = System.getProperty("user.name");
        JFrame frame = new JFrame(username);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(new MyCoreInputPanel());
        frame.setSize(new Dimension(600,700));
        frame.setVisible(true);
    }
}
