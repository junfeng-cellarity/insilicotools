package com.insilico.application.insilicotools.gui;


import chemaxon.calculations.clean.Cleaner;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.InHouseCollectionDAO;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.sql.SQLException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.regex.Pattern;

public class PasteSmilesDialog extends StandardDialog{
    private JTextArea textArea;
    Vector<PropertyMolecule> molecules;
    JPanel mainPanel;
    ButtonPanel buttonPanel;
    JProgressBar progressBar;
    boolean generateChiral = false;


    public PasteSmilesDialog() {
        super();
        setModal(true);
        setSize(new Dimension(800,600));
        initializeComponents();
    }

    public Vector<PropertyMolecule> getMolecules() {
        return molecules;
    }

    protected void initializeComponents() {
        progressBar = new JProgressBar(0,100);
        molecules = new Vector<PropertyMolecule>();
        mainPanel = new JPanel(new BorderLayout());
        mainPanel.setBorder(new TitledBorder("Paste Smiles/CY_NUMBER Here"));
        textArea = new JTextArea();
        mainPanel.add(new JScrollPane(textArea), BorderLayout.CENTER);
        mainPanel.add(progressBar, BorderLayout.SOUTH);

        buttonPanel = new ButtonPanel();
        buttonPanel.setAlignment(SwingConstants.CENTER);
        final JButton clearButton = new JButton("Clear", new ImageIcon(Toolkit.getDefaultToolkit().createImage(getClass().getResource("/toolbarButtonGraphics/general/Delete16.gif"))));
        clearButton.setEnabled(false);
        clearButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                textArea.setText("");
            }
        });
        buttonPanel.add(clearButton);

        final JButton cancelButton = new JButton("Cancel");
        cancelButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                setDialogResult(StandardDialog.RESULT_CANCELLED);
                setVisible(false);
            }
        });
        buttonPanel.add(cancelButton);

        final JButton loadSmilesButton = new JButton("Load Smiles/CY", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Import24.gif")));
        loadSmilesButton.setEnabled(false);
        loadSmilesButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                final String text = textArea.getText();
                progressBar.setValue(1);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected void process(List chunks) {
                        Vector v = (Vector) chunks.get(chunks.size() - 1);
                        String note = (String) v.get(0);
                        int progress = (Integer) v.get(1);
                        progressBar.setValue(progress);
                        progressBar.setString(note);
                    }

                    @Override
                    protected Object doInBackground() throws Exception {
                        return parseCYNumberorSmilesToDesignMolecules(text, generateChiral, new ProgressReporter() {
                            @Override
                            public void reportProgress(String note, int progress) {
                                Vector v = new Vector();
                                v.add(note);
                                v.add(progress);
                                publish(v);
                            }
                        });
                    }

                    @Override
                    protected void done() {
                        try {
                            ArrayList molList = (ArrayList) get();
                            if (molList.size() > 0) {
                                molecules.clear();
                                for (Object m : molList) {
                                    molecules.add((PropertyMolecule)m);
                                }
                                //ChemFunc.calculateOEProperty(molecules);
                                firePropertyChange("moleculeLoaded",true,false);
                            } else {
                                JOptionPane.showMessageDialog(PasteSmilesDialog.this, "No molecule found.");
                            }
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(PasteSmilesDialog.this, e1.getMessage());
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(PasteSmilesDialog.this, e1.getMessage());
                        }finally {
                            progressBar.setString("");
                            progressBar.setValue(0);
                        }
                    }
                };
                sw.addPropertyChangeListener(new PropertyChangeListener() {
                    @Override
                    public void propertyChange(PropertyChangeEvent evt) {
                        if(evt.getPropertyName().equals("moleculeLoaded")){
                            setDialogResult(StandardDialog.RESULT_AFFIRMED);
                            setVisible(false);
                        }
                    }
                });
                sw.execute();
            }
        });
        buttonPanel.add(loadSmilesButton);

        textArea.getDocument().addDocumentListener(new DocumentListener() {
            public void changedUpdate(DocumentEvent e) {
                boolean b = textArea.getText().trim().length() > 0;
                loadSmilesButton.setEnabled(b);
                clearButton.setEnabled(b);
            }

            public void insertUpdate(DocumentEvent e) {
                boolean b = textArea.getText().trim().length() > 0;
                loadSmilesButton.setEnabled(b);
                clearButton.setEnabled(b);
            }

            public void removeUpdate(DocumentEvent e) {
                boolean b = textArea.getText().trim().length() > 0;
                loadSmilesButton.setEnabled(b);
                clearButton.setEnabled(b);
            }
        });

        final JPopupMenu popup = new JPopupMenu();

        JMenuItem cutMenuItem = new JMenuItem("Cut", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Cut16.gif")));
        cutMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                textArea.cut();
            }
        });
        popup.add(cutMenuItem);

        JMenuItem copyMenuItem = new JMenuItem("Copy", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Copy16.gif")));
        copyMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                textArea.copy();
            }
        });
        popup.add(copyMenuItem);

        JMenuItem pasteMenuItem = new JMenuItem("Paste", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Paste16.gif")));
        pasteMenuItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                textArea.paste();
            }
        });
        popup.add(pasteMenuItem);

        final JCheckBoxMenuItem enumeratItem = new JCheckBoxMenuItem("Enumerate all stereoisomers",false);
        enumeratItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                generateChiral = enumeratItem.isSelected();
            }
        });
        popup.add(enumeratItem);

        textArea.addMouseListener(new MouseAdapter() {
            public void mousePressed(MouseEvent e) {
                if (e.isPopupTrigger()) popup.show(e.getComponent(), e.getX(), e.getY());
            }

            public void mouseReleased(MouseEvent e) {
                if (e.isPopupTrigger()) popup.show(e.getComponent(), e.getX(), e.getY());
            }
        });
    }

    @Override
    public JComponent createBannerPanel() {
        return null;
    }

    @Override
    public JComponent createContentPanel() {
        return mainPanel;
    }

    @Override
    public ButtonPanel createButtonPanel() {
        return buttonPanel;
    }

    public void setVisible(boolean b) {
        super.setVisible(b);

        if (b) {
            textArea.requestFocus();
        }
    }


    private static String convertWrongLotNumberToRightLotNumber(String wrongLotNumber){
        String[] s = wrongLotNumber.split("-", 3);
        return String.format("%s-%s-0%s",s[0],s[1],s[2]);
        //add 0 to lot number
        //assuming the
    }

    private static ArrayList<PropertyMolecule> parseCYNumberorSmilesToDesignMolecules(String text, boolean generateChiral, ProgressReporter progressReporter) {
        Pattern cynumberRe = Pattern.compile("^CY-[0-9][0-9][0-9][0-9][0-9][0-9][0-9]$");
        Pattern cynumberReWithLot = Pattern.compile("^CY-[0-9][0-9][0-9][0-9][0-9][0-9][0-9]-[0-9][0-9]$");

        ArrayList<PropertyMolecule> molecules = new ArrayList<PropertyMolecule>();
        String[] lines = text.split("\\n");
        if(lines.length>=200) {
            Vector<String> bioList = new Vector<String>();
            Vector<String> nameList = new Vector<String>();
            for (String line : lines) {
                line = line.trim();
                if (cynumberRe.matcher(line).matches()) {
                    bioList.add(line);
                }else if(cynumberReWithLot.matcher(line).matches()){
                    String cy_number = line.substring(0,line.lastIndexOf("-"));
                    bioList.add(cy_number);
                }
                else {
                    nameList.add(line);
                }
            }
            try {
                final HashMap<String,Integer> idxDict = new HashMap<String, Integer>();
                int idx = 0;
                for(String bio:bioList){
                    idxDict.put(bio,idx++);
                }
                Vector<PropertyMolecule> cmpds = InHouseCollectionDAO.getInstance().getMolFromCYNumberBatch(bioList, progressReporter);
                Collections.sort(cmpds, new Comparator<PropertyMolecule>() {
                    @Override
                    public int compare(PropertyMolecule o1, PropertyMolecule o2) {
                        return idxDict.get(o1.getName()).compareTo(idxDict.get(o2.getName()));
                    }
                });
                molecules.addAll(cmpds);
            } catch (SQLException e) {
                e.printStackTrace();
            }

            for (int i = 0; i < nameList.size(); i++) {
                progressReporter.reportProgress("Loading molecule from smiles ...", 100 * i / nameList.size());
                String line = nameList.get(i).trim().replaceAll("\"", "").replaceAll("\'", "");
                if (line.length() > 0) {
                    String[] args = line.split("\\s+");
                    OEGraphMol mol = new OEGraphMol();
                    boolean success = oechem.OEParseSmiles(mol, args[0]);
                    if (!success) {
                        System.err.println(line + " not imported.");
                        continue;
                    }

                    if (mol.NumAtoms() == 0) {
                        System.err.println(args[0] + " not imported.");
                        continue;
                    }
                    if (args.length >= 2) {
                        mol.SetTitle(args[1]);
                    }
                    processPropertyMolecules(generateChiral,molecules,mol);
                }
            }
            return molecules;
        }

        for (int i = 0; i < lines.length; i++) {
            progressReporter.reportProgress("Loading molecules ...",100*i/lines.length);
            String line = lines[i].trim().replaceAll("\"","").replaceAll("\'","");
            if (line.length() > 0) {
                String[] args = line.split("\\s+");
                String cy_number = args[0].toUpperCase();
                if(cynumberRe.matcher(cy_number).matches()){
                    try {
                        OEGraphMol mol = InHouseCollectionDAO.getInstance().getMolFromCYNumber(cy_number);
                        if(mol!=null&&mol.NumAtoms()>0){
                            molecules.add(new PropertyMolecule(mol));
                        }
                    } catch (SQLException e) {
                        e.printStackTrace();
                    }
                }else if(cynumberReWithLot.matcher(line).matches()){
                    cy_number = line.substring(0,line.lastIndexOf("-"));
                    System.out.println(cy_number);
                    try {
                        OEGraphMol mol = InHouseCollectionDAO.getInstance().getMolFromCYNumber(cy_number);
                        if(mol!=null&&mol.NumAtoms()>0){
                            molecules.add(new PropertyMolecule(mol));
                        }
                    } catch (SQLException e) {
                        e.printStackTrace();
                    }
                }
                else if(cy_number.startsWith("CY")){
                    cy_number = ChemFunc.formatCYNumber(cy_number);
                    try {
                        OEGraphMol mol = InHouseCollectionDAO.getInstance().getMolFromCYNumber(cy_number);
                        if(mol!=null&&mol.NumAtoms()>0){
                            molecules.add(new PropertyMolecule(mol));
                        }
                    } catch (SQLException e) {
                        e.printStackTrace();
                    }
                }
                else{
                    OEGraphMol mol = new OEGraphMol();
                    boolean success = oechem.OEParseSmiles(mol, args[0]);
                    if(!success){
                        System.err.println(line+ " not imported.");
                        continue;
                    }

                    if(mol.NumAtoms()==0){
                        System.err.println(args[0]+" not imported.");
                        continue;
                    }
                    if (args.length >= 2) {
                        mol.SetTitle(args[1]);
                    }
                    processPropertyMolecules(generateChiral, molecules, mol);
                }
            }
        }
        return molecules;
    }

    public static void processPropertyMolecules(boolean generateChiral, List<PropertyMolecule> molecules, OEGraphMol mol) {
        if(generateChiral){
            Vector<OEGraphMol> chiralMols = OEChemFunc.getChiralMols(mol, true);
            for(OEGraphMol m:chiralMols){
                m.SetDimension(1);
                System.out.println(oechem.OEMolToSmiles(m));
                Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(m);
                Cleaner.clean(molecule,2);
                m = OEChemFunc.getInstance().convertChemAxonMol(molecule);
                System.out.println(oechem.OEMolToSmiles(m));
                molecules.add(new PropertyMolecule(m));
            }
        }else {
            if(mol.GetDimension()!=2){
                Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(mol);
                Cleaner.clean(molecule,2);
                mol = OEChemFunc.getInstance().convertChemAxonMol(molecule);
            }
            molecules.add(new PropertyMolecule(mol));
        }
    }

    public static void main(String[] args) {
        PasteSmilesDialog dialog = new PasteSmilesDialog();
        dialog.setVisible(true);
    }

}