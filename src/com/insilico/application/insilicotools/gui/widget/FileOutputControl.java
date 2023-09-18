package com.insilico.application.insilicotools.gui.widget;

import com.jidesoft.swing.PartialEtchedBorder;

import javax.swing.*;
import javax.swing.filechooser.FileFilter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

/**
 * Created by jfeng1 on 2/28/17.
 */
public class FileOutputControl extends JPanel {
    File file;
    JTextField fileField;
    JFileChooser fc;
    File currentDirectory;

    public String getPropertyName(){
        return "FileLoaded";
    }


    public File getCurrentDirectory() {
        return currentDirectory;
    }

    public void setCurrentDirectory(File currentDirectory) {
        this.currentDirectory = currentDirectory;
    }

    public FileOutputControl(final File currentDirectory){
        this(currentDirectory,".sdf");
    }

    public FileOutputControl(final File currentDirectory1, final String suffix) {
        super(new FlowLayout(FlowLayout.LEADING));
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p,BoxLayout.LINE_AXIS));
        setPreferredSize(new Dimension(600,45));
        fileField = new JTextField();
        fileField.setEnabled(false);
        fileField.setPreferredSize(new Dimension(500,35));
        this.currentDirectory = currentDirectory1;
        fc = new JFileChooser();
        fc.setFileFilter(new FileFilter() {
            @Override
            public boolean accept(File f) {
                if(f.isDirectory()){
                    return true;
                }
                if(f.getName().endsWith(suffix)){
                    return true;
                }
                return false;
            }

            @Override
            public String getDescription() {
                return suffix;
            }
        });
        JButton btn = new JButton("...");
        btn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(currentDirectory!=null){
                    fc.setCurrentDirectory(currentDirectory);
                }
                int option = fc.showSaveDialog(FileOutputControl.this);
                if(option==JFileChooser.APPROVE_OPTION){
                    file = fc.getSelectedFile();
                    if(!file.getName().endsWith(suffix)){
                        file = new File(file.getAbsolutePath()+suffix);
                    }
                    if(file.exists()){
                        int reply = JOptionPane.showConfirmDialog(FileOutputControl.this,"File exists, Overwrite?");
                        if(reply==JOptionPane.NO_OPTION){
                            return;
                        }
                    }
                    fileField.setText(file.getAbsolutePath());
                    if(file.getParentFile()!=null&&file.getParentFile().isDirectory()) {
                        currentDirectory = file.getParentFile();
                    }
                    FileOutputControl.this.firePropertyChange(getPropertyName(),true,false);
                }
            }
        });
        p.add(btn);
        p.add(fileField);
        add(p);
        setBorder(new PartialEtchedBorder());
    }

    public File getFile() {
        return file;
    }


    public static void main(String[] args) {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(new FileOutputControl(new File(System.getProperty("user.dir"))));
        frame.pack();
        frame.setVisible(true);
    }
}
