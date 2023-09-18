package com.insilico.application.insilicotools.cmdline;

import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.PropertySelectionPanel;
import com.insilico.application.insilicotools.gui.widget.FileInputControl;
import com.insilico.application.insilicotools.gui.widget.FileOutputControl;
import com.insilico.application.insilicotools.util.ChemFunc;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutionException;

public class BatchDescriptorGenerationPanel extends JPanel {
    File currentDirectory = new File(System.getProperty("user.dir"));
    FileInputControl input_file_control;
    FileOutputControl output_file_control;
    PropertySelectionPanel mainPanel;
    File input_file;
    File output_file;
    boolean isSubmitted = false;
    Vector<String> properties = new Vector<>();

    JButton okBtn = new JButton("OK");
    JButton cancelBtn = new JButton("Cancel");
    JProgressBar progressBar = new JProgressBar(0,100);

    public BatchDescriptorGenerationPanel() {
        super(new BorderLayout());
        okBtn.setEnabled(false);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(properties.isEmpty()){
                    JOptionPane.showMessageDialog(BatchDescriptorGenerationPanel.this,"No properties selected!");
                    return;
                }
                isSubmitted = true;
                SwingWorker sw = new SwingWorker() {
                    @Override
                    protected void process(List chunks) {
                        int progress = (Integer)chunks.get(chunks.size() - 1);
                        if(progress==-1){
                            progressBar.setIndeterminate(true);
                        }else {
                            progressBar.setValue(progress);
                        }
                    }

                    @Override
                    protected Object doInBackground() throws Exception {
                        BatchJobs.calculateProperties_all(input_file.getAbsolutePath(), output_file.getAbsolutePath(), properties, new ProgressReporter() {
                            @Override
                            public void reportProgress(String note, final int progress) {
                                publish(progress);
                            }
                        });
                        return null;
                    }

                    @Override
                    protected void done() {
                        try {
                            get();
                            BatchDescriptorGenerationPanel.this.getTopLevelAncestor().setVisible(false);
                        } catch (InterruptedException | ExecutionException e1) {
                            e1.printStackTrace();
                            JOptionPane.showMessageDialog(BatchDescriptorGenerationPanel.this,e1.getMessage());
                        }
                    }
                };
                sw.execute();
            }
        });
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isSubmitted = false;
                BatchDescriptorGenerationPanel.this.getTopLevelAncestor().setVisible(false);
            }
        });
        JPanel gridPanel = new JPanel(new GridLayout(2,1));
        input_file_control = new FileInputControl(currentDirectory);
        output_file_control = new FileOutputControl(currentDirectory);
        JPanel inputPanel = new JPanel();
        inputPanel.add(new JLabel("Input File:"));
        inputPanel.add(input_file_control);
        input_file_control.addPropertyChangeListener(input_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                input_file = input_file_control.getFile();
                output_file_control.setCurrentDirectory(input_file_control.getCurrentDirectory());
                properties.clear();
                properties.addAll(mainPanel.getSelectedProperties());
                checkCondition();
            }
        });

        gridPanel.add(inputPanel);
        JPanel outputPanel = new JPanel();
        outputPanel.add(new JLabel("Output File:"));
        outputPanel.add(output_file_control);
        output_file_control.addPropertyChangeListener(output_file_control.getPropertyName(), new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                output_file = output_file_control.getFile();
                properties.clear();
                properties.addAll(mainPanel.getSelectedProperties());
                checkCondition();
            }
        });

        gridPanel.add(outputPanel);
        mainPanel = new PropertySelectionPanel(ChemFunc.properties,ChemFunc.adme_models,ChemFunc.safety_models,ChemFunc.other_models);
        mainPanel.addPropertyChangeListener("selection", new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                properties.clear();
                properties.addAll(mainPanel.getSelectedProperties());
                checkCondition();
            }
        });
        add(gridPanel,BorderLayout.NORTH);
        add(mainPanel,BorderLayout.CENTER);
        JPanel gridPanel2 = new JPanel(new GridLayout(2,1));
        JPanel btnPanel = new JPanel();
        btnPanel.add(okBtn);
        btnPanel.add(cancelBtn);
        gridPanel2.add(btnPanel);
        gridPanel2.add(progressBar);
        add(gridPanel2,BorderLayout.SOUTH);

    }

    private void checkCondition(){
        if(input_file!=null&&output_file!=null&&properties.size()>0){
            okBtn.setEnabled(true);
        }else{
            okBtn.setEnabled(false);
        }
    }

    public static void main(String[] args) {
        JDialog dialog = new JDialog();
        BatchDescriptorGenerationPanel p = new BatchDescriptorGenerationPanel();
        dialog.getContentPane().add(p);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setSize(new Dimension(800,600));
        dialog.setVisible(true);
    }

}
