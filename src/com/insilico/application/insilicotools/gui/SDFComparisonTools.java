package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.gui.widget.FileInputControl;
import com.insilico.application.insilicotools.gui.widget.FileOutputControl;
import openeye.oechem.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 4/27/17.
 */
public class SDFComparisonTools extends JPanel {
    File currentDirectory = new File(System.getProperty("user.dir"));
    FileInputControl origin_file_control;
    FileInputControl filter_file_control;
    FileOutputControl uniq_file_control;
    FileOutputControl duplicate_file_control;
    File origin_file;
    File filter_file;
    File uniq_file;
    File duplicate_file;
    boolean isSubmitted = false;

    JButton okBtn;
    JButton cancelBtn;
    JProgressBar progressBar = new JProgressBar(0,100);

    private void checkCondition(){
        if(origin_file!=null&&filter_file!=null&&uniq_file!=null&&duplicate_file!=null){
            okBtn.setEnabled(true);
        }else{
            okBtn.setEnabled(false);
        }

    }

    public SDFComparisonTools() {
        super(new BorderLayout());

        JPanel p = new JPanel();
        BoxLayout layout = new BoxLayout(p, BoxLayout.Y_AXIS);
        p.setLayout(layout);
        origin_file_control = new FileInputControl(currentDirectory);
        origin_file_control.addPropertyChangeListener(origin_file_control.getPropertyName(),new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                currentDirectory = origin_file_control.getCurrentDirectory();
                filter_file_control.setCurrentDirectory(currentDirectory);
                uniq_file_control.setCurrentDirectory(currentDirectory);
                duplicate_file_control.setCurrentDirectory(currentDirectory);
                origin_file = origin_file_control.getFile();
                checkCondition();
            }
        });
        filter_file_control = new FileInputControl(currentDirectory);
        filter_file_control.addPropertyChangeListener(filter_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                filter_file = filter_file_control.getFile();
                checkCondition();
            }
        });
        uniq_file_control = new FileOutputControl(currentDirectory);
        uniq_file_control.addPropertyChangeListener(uniq_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                uniq_file = uniq_file_control.getFile();
                checkCondition();
            }
        });
        duplicate_file_control = new FileOutputControl(currentDirectory);
        duplicate_file_control.addPropertyChangeListener(duplicate_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                duplicate_file = duplicate_file_control.getFile();
                checkCondition();
            }
        });
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Input SDF 1:"));
        p1.add(origin_file_control);
        JPanel p2 = new JPanel();
        p2.add(new JLabel("Input SDF 2:"));
        p2.add(filter_file_control);
        JPanel p3 = new JPanel();
        p3.add(new JLabel("Unique SDF:"));
        p3.add(uniq_file_control);
        JPanel p4 = new JPanel();
        p4.add(new JLabel("Duplicate SDF:"));
        p4.add(duplicate_file_control);
        p.add(p1);
        p.add(p2);
        p.add(p3);
        p.add(p4);
        progressBar.setBorderPainted(true);
        progressBar.setBorder(new TitledBorder("Progress"));
        p.add(progressBar);

        add(p,BorderLayout.CENTER);

        JPanel btnPanel = new JPanel();
        okBtn = new JButton("OK");
        okBtn.setEnabled(false);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(origin_file==null||!origin_file.canRead()||filter_file==null||!filter_file.canRead()){
                    JOptionPane.showMessageDialog(SDFComparisonTools.this,"Input file is not readable!");
                    return;
                }
                if(uniq_file==null|duplicate_file==null){
                    JOptionPane.showMessageDialog(SDFComparisonTools.this,"No output file is specified!");
                    return;
                }
                isSubmitted = true;
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        HashMap<String,Integer> map = new HashMap<String, Integer>();
                        OEGraphMol mol = new OEGraphMol();
                        oemolistream ifs_0 = new oemolistream(origin_file.getAbsolutePath());
                        oemolistream ifs_1 = new oemolistream(origin_file.getAbsolutePath());
                        OEMolDatabase mdb = new OEMolDatabase(ifs_0);
                        ifs_0.close();
                        int numMols = mdb.GetMaxMolIdx();
                        oemolistream ifs_2 = new oemolistream(filter_file.getAbsolutePath());
                        oemolostream ofs_1 = new oemolostream(uniq_file.getAbsolutePath());
                        oemolostream ofs_2 = new oemolostream(duplicate_file.getAbsolutePath());
                        publish(1);
                        while(oechem.OEReadMolecule(ifs_2,mol)){
                            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                            map.put(oechem.OEMolToSmiles(mol),1);
                        }
                        ifs_2.close();
                        int n = 0;
                        while(oechem.OEReadMolecule(ifs_1,mol)){
                            n += 1;
                            publish(100*n/numMols);
                            OEGraphMol mol2 = new OEGraphMol(mol);
                            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol2);
                            String smiles = oechem.OEMolToSmiles(mol2);
                            if(map.containsKey(smiles)){
                                oechem.OEWriteMolecule(ofs_2,mol);
                            }else{
                                oechem.OEWriteMolecule(ofs_1,mol);
                            }
                        }
                        ifs_1.close();
                        ofs_1.close();
                        ofs_2.close();
                        return null;
                    }

                    @Override
                    protected void process(List chunks) {
                        int progress = (Integer)chunks.get(chunks.size()-1);
                        progressBar.setValue(progress);
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            SDFComparisonTools.this.getTopLevelAncestor().setVisible(false);
                            JOptionPane.showMessageDialog(SDFComparisonTools.this,"Finished.");
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                        }
                    }
                };
                okBtn.setEnabled(false);
                sw.execute();
            }
        });
        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isSubmitted = false;
                SDFComparisonTools.this.getTopLevelAncestor().setVisible(false);
            }
        });
        btnPanel.add(okBtn);
        btnPanel.add(cancelBtn);
        add(btnPanel,BorderLayout.SOUTH);
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(new SDFComparisonTools());
        f.pack();
        f.setVisible(true);
    }
}
