package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.data.Chemist;
import com.insilico.application.insilicotools.data.Project;
import com.insilico.application.insilicotools.data.Status;
import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.gui.widget.FileInputControl;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import com.jidesoft.swing.PartialEtchedBorder;
import openeye.oechem.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.Vector;

/**
 * Created by jfeng1 on 2/9/17.
 */
public class MyCompoundsBatchRegistrationDialog extends StandardDialog {
    JPanel mainPanel;
    FileInputControl fileInputControl;
    Vector<String> nameTags;
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;
    Status MyStatus = null;
    Project MyProject = null;
    String nameTag = null;

    String sdfFilePath;
    JComboBox projectCB;
    JComboBox<Chemist> chemistCB = new JComboBox<>();
    int chemist_id = 0;

    public String getSdfFilePath() {
        return sdfFilePath;
    }

    public MyCompoundsBatchRegistrationDialog() {
        super();
        setModal(true);
        setPreferredSize(new Dimension(750,250));
        initializeComponents();
        pack();
    }

    public void setMyProject(Project MyProject) {
        if(MyProject!=null) {
            this.MyProject = MyProject;
            projectCB.setSelectedItem(MyProject);
        }
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public void initializeComponents() {
        mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel,BoxLayout.PAGE_AXIS));
        nameTags = new Vector<String>();
        File currentDirectory = new File(System.getProperty("user.dir"));
        fileInputControl = new FileInputControl(currentDirectory);

        JPanel p1 = new JPanel(new FlowLayout(FlowLayout.LEADING));
        p1.add(new JLabel("1: Load SDF file:"));
        p1.add(fileInputControl);

        JPanel p2 = new JPanel(new FlowLayout(FlowLayout.LEADING));
        p2.add(new JLabel("2:Select name tag:"));
        final JComboBox nameTagCb = new JComboBox();
        nameTagCb.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                nameTag = (String) nameTagCb.getSelectedItem();
            }
        });
        p2.add(nameTagCb);

        final Vector<Project> My_projects = FrontierDAO.getInstance().getMy_projects();
        MyProject = My_projects.isEmpty()?null:My_projects.get(0);
        projectCB = new JComboBox(My_projects);
        projectCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MyProject = (Project) projectCB.getSelectedItem();
                Vector<Chemist> chemistsByProject = FrontierDAO.getInstance().getChemistsByProject(MyProject.getProject_id());
                if(chemistsByProject!=null) {
                    chemistCB.setModel(new DefaultComboBoxModel(chemistsByProject));
                }else{
                    chemistCB.setModel(new DefaultComboBoxModel());
                }
            }
        });
        Vector<Status> My_status = FrontierDAO.getInstance().getMy_status();
        MyStatus = My_status.isEmpty()?null:My_status.get(0);
        final JComboBox statusCB = new JComboBox(My_status);
        statusCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                MyStatus = (Status) statusCB.getSelectedItem();
            }
        });
        JPanel p3 = new JPanel(new FlowLayout(FlowLayout.LEADING));
        p3.add(new JLabel("3:Project:"));
        p3.add(projectCB);
        p3.add(new JLabel("Status:"));
        p3.add(statusCB);
        p3.add(new JLabel("Assign To:"));
        p3.add(chemistCB);

        fileInputControl.addPropertyChangeListener("FileLoaded",new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                File file = fileInputControl.getFile();
                if(file !=null&&file.canRead()){
                    String sdf = file.getAbsolutePath();
                    oemolistream ifs = new oemolistream();
                    ifs.open(sdf);
                    nameTags.clear();
                    OEGraphMol mol = new OEGraphMol();
                    int max = 10;
                    int n = 0;
                    while(oechem.OEReadMolecule(ifs,mol)&&n<max){
                        OESDDataIter pairs = oechem.OEGetSDDataPairs(mol);
                        while(pairs.iterator().hasNext()){
                            OESDDataPair next = pairs.iterator().next();
                            String tag = next.GetTag();
                            if(!nameTags.contains(tag)){
                                nameTags.add(tag);
                            }
                            n += 1;
                        }
                    }
                    if(!nameTags.isEmpty()){
                        nameTags.insertElementAt("",0);
                        nameTag = nameTags.get(0);
                        nameTagCb.setEnabled(true);
                        nameTagCb.setModel(new DefaultComboBoxModel(nameTags));
                        okBtn.setEnabled(true);
                    }else{
                        nameTag = null;
                        nameTagCb.setEnabled(false);
                        nameTagCb.setModel(new DefaultComboBoxModel());
                        okBtn.setEnabled(true);
                    }

                }
            }
        });

        mainPanel.add(p1);
        mainPanel.add(p2);
        mainPanel.add(p3);
        mainPanel.setBorder(new PartialEtchedBorder());

        okBtn = new JButton("OK");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                sdfFilePath = fileInputControl.getFile().getAbsolutePath();
                isCommitted = true;
                setVisible(false);
            }
        });
        okBtn.setEnabled(false);
        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                setVisible(false);
            }
        });
    }


    public static void main(String[] args) {
        MyCompoundsBatchRegistrationDialog dialog = new MyCompoundsBatchRegistrationDialog();
        dialog.setVisible(true);

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
        ButtonPanel p = new ButtonPanel(SwingConstants.CENTER);
        p.add(okBtn,ButtonPanel.AFFIRMATIVE_BUTTON);
        p.add(cancelBtn,ButtonPanel.CANCEL_BUTTON);
        return p;
    }

    public String getMolNameTag(){
        return nameTag;
    }

    public Status getMyStatus() {
        return MyStatus;
    }

    public Project getMyProject() {
        return MyProject;
    }

    public Chemist getMyChemist(){
        if(chemistCB.getModel().getSize()==0){
            return null;
        }else{
            return chemistCB.getItemAt(chemistCB.getSelectedIndex());
        }
    }
}
