package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.gui.widget.FileInputControl;
import com.insilico.application.insilicotools.gui.widget.FileOutputControl;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.net.URI;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutionException;

/**
 * Created by jfeng1 on 4/27/17.
 */
public class NMRConversionTool extends JPanel {
    File currentDirectory = new File(System.getProperty("user.dir"));
    FileInputControl input_file_control;
    FileOutputControl output_file_control;
    File input_file;
    File output_file;
    boolean isSubmitted = false;

    JButton okBtn;
    JButton cancelBtn;
    JProgressBar progressBar = new JProgressBar(0,100);

    private void checkCondition(){
        if(input_file!=null&&output_file!=null){
            okBtn.setEnabled(true);
        }else{
            okBtn.setEnabled(false);
        }

    }

    public NMRConversionTool() {
        super(new BorderLayout());

        JPanel p = new JPanel();
        BoxLayout layout = new BoxLayout(p, BoxLayout.Y_AXIS);
        p.setLayout(layout);
        input_file_control = new FileInputControl(currentDirectory,".pdb");
        input_file_control.addPropertyChangeListener(input_file_control.getPropertyName(),new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                currentDirectory = input_file_control.getCurrentDirectory();
                output_file_control.setCurrentDirectory(currentDirectory);
                input_file = input_file_control.getFile();
                checkCondition();
            }
        });
        output_file_control = new FileOutputControl(currentDirectory,".lib");
        output_file_control.addPropertyChangeListener(output_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                output_file = output_file_control.getFile();
                checkCondition();
            }
        });
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Input PDB:"));
        p1.add(input_file_control);
        JPanel p3 = new JPanel();
        p3.add(new JLabel("Output Lib:"));
        p3.add(output_file_control);
        p.add(p1);
        p.add(p3);
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
                if(input_file==null||!input_file.canRead()){
                    JOptionPane.showMessageDialog(NMRConversionTool.this,"Input file is not readable!");
                    return;
                }
                if(output_file==null){
                    JOptionPane.showMessageDialog(NMRConversionTool.this,"No output file is specified!");
                    return;
                }
                isSubmitted = true;
                progressBar.setIndeterminate(true);
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected Object doInBackground() throws Exception {
                        String pdbString = FileUtil.readFileToString(input_file);
                        String libString = ChemFunc.convertNMRLib(pdbString);
                        FileUtil.writeStringToFile(output_file, libString);
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
                            NMRConversionTool.this.getTopLevelAncestor().setVisible(false);
                            JOptionPane.showMessageDialog(NMRConversionTool.this,"Finished.");
                        } catch (InterruptedException e1) {
                            e1.printStackTrace();
                        } catch (ExecutionException e1) {
                            e1.printStackTrace();
                        }finally {
                            progressBar.setIndeterminate(false);
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
                NMRConversionTool.this.getTopLevelAncestor().setVisible(false);
            }
        });
        btnPanel.add(okBtn);
        btnPanel.add(cancelBtn);
        add(btnPanel,BorderLayout.SOUTH);
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(new NMRConversionTool());
        f.pack();
        f.setVisible(true);
    }
}
