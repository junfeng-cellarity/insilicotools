package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.gui.widget.FileInputControl;
import com.insilico.application.insilicotools.gui.widget.FileOutputControl;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import openeye.oechem.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 4/27/17.
 */
public class ChiralEnumerationTools extends JPanel {
    File currentDirectory = new File(System.getProperty("user.dir"));
    FileInputControl origin_file_control;
    FileOutputControl chiral_file_control;
    File origin_file;
    File chiral_file;
    boolean isSubmitted = false;

    JButton okBtn;
    JButton cancelBtn;
    JProgressBar progressBar = new JProgressBar(0,100);

    private void checkCondition(){
        if(origin_file!=null&&chiral_file!=null){
            okBtn.setEnabled(true);
        }else{
            okBtn.setEnabled(false);
        }

    }

    public ChiralEnumerationTools() {
        super(new BorderLayout());

        JPanel p = new JPanel();
        BoxLayout layout = new BoxLayout(p, BoxLayout.Y_AXIS);
        p.setLayout(layout);
        origin_file_control = new FileInputControl(currentDirectory);
        origin_file_control.addPropertyChangeListener(origin_file_control.getPropertyName(),new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                currentDirectory = origin_file_control.getCurrentDirectory();
                chiral_file_control.setCurrentDirectory(currentDirectory);
                origin_file = origin_file_control.getFile();
                checkCondition();
            }
        });
        chiral_file_control = new FileOutputControl(currentDirectory);
        chiral_file_control.addPropertyChangeListener(chiral_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                chiral_file = chiral_file_control.getFile();
                checkCondition();
            }
        });
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Input SDF:"));
        p1.add(origin_file_control);
        JPanel p2 = new JPanel();
        p2.add(new JLabel("Chiral Output SDF:"));
        p2.add(chiral_file_control);
        p.add(p1);
        p.add(p2);
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
                if(origin_file==null||!origin_file.canRead()){
                    JOptionPane.showMessageDialog(ChiralEnumerationTools.this,"Input file is not readable!");
                    return;
                }
                if(chiral_file==null){
                    JOptionPane.showMessageDialog(ChiralEnumerationTools.this,"No output file is specified!");
                    return;
                }
                isSubmitted = true;
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        OEGraphMol mol = new OEGraphMol();
                        oemolistream tmp_ifs = new oemolistream(origin_file.getAbsolutePath());
                        OEMolDatabase mdb = new OEMolDatabase(tmp_ifs);
                        tmp_ifs.close();
                        int numMols = mdb.GetMaxMolIdx();
                        oemolistream ifs = new oemolistream(origin_file.getAbsolutePath());
                        oemolostream ofs = new oemolostream(chiral_file.getAbsolutePath());
                        publish(1);
                        while(oechem.OEReadMolecule(ifs,mol)){
                            Vector<OEGraphMol> chiralMols = OEChemFunc.getChiralMols(mol, true);
                            int n =0;
                            for(OEGraphMol chiralMol:chiralMols){
                                oechem.OECopySDData(mol,chiralMol);
                                chiralMol.SetTitle(String.format("%s_E%d",mol.GetTitle(),++n));
                                System.out.println(oechem.OEMolToSmiles(chiralMol));
                            }
                            for(int i=0;i<4;i++){
                                if(i<chiralMols.size()){
                                    OEGraphMol mol1 = chiralMols.get(i);
                                    mol1 = OEChemFunc.getInstance().getMol2D(mol1);
                                    oechem.OEWriteMolecule(ofs, mol1);
                                }else{
                                    OEGraphMol tmpMol = new OEGraphMol();
                                    oechem.OECopySDData(mol,tmpMol);
                                    oechem.OEWriteMolecule(ofs,tmpMol);
                                    System.out.println(ChemFunc.getMolString(tmpMol));
                                }
                            }

                        }
                        ifs.close();
                        ofs.close();
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
                            ChiralEnumerationTools.this.getTopLevelAncestor().setVisible(false);
                            JOptionPane.showMessageDialog(ChiralEnumerationTools.this,"Finished.");
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
                ChiralEnumerationTools.this.getTopLevelAncestor().setVisible(false);
            }
        });
        btnPanel.add(okBtn);
        btnPanel.add(cancelBtn);
        add(btnPanel,BorderLayout.SOUTH);
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(new ChiralEnumerationTools());
        f.pack();
        f.setVisible(true);
    }
}
