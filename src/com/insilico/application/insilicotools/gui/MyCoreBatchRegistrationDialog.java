package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.gui.widget.FileInputControl;
import com.google.common.base.Strings;
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
public class MyCoreBatchRegistrationDialog extends StandardDialog {
    JPanel mainPanel;
    FileInputControl fileInputControl;
    Vector<String> nameTags;
    JButton okBtn;
    JButton cancelBtn;
    boolean isCommitted = false;

    String libraryNameTag;
    String sdfFilePath;
    String source = "Cellarity";

    public String getLibraryNameTag() {
        return libraryNameTag;
    }

    public String getSdfFilePath() {
        return sdfFilePath;
    }

    public String getSource() {
        return source;
    }

    public MyCoreBatchRegistrationDialog() {
        super();
        setModal(true);
        setSize(new Dimension(650,250));
        initializeComponents();
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
        nameTagCb.setEnabled(false);
        p2.add(nameTagCb);

        JPanel p3 = new JPanel(new FlowLayout(FlowLayout.LEADING));
        p3.add(new JLabel("3:Add source: (e.g. My)"));
        final JTextField textField = new JTextField(20);
        textField.setText("My");
        textField.setEnabled(false);
        p3.add(textField);

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
                        nameTagCb.setEnabled(true);
                        nameTagCb.setModel(new DefaultComboBoxModel(nameTags));
                        textField.setEnabled(true);
                        okBtn.setEnabled(true);
                    }else{
                        nameTagCb.setEnabled(false);
                        nameTagCb.setModel(new DefaultComboBoxModel());
                        okBtn.setEnabled(false);
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
                libraryNameTag = (String)nameTagCb.getSelectedItem();
                source = textField.getText();
                if(Strings.isNullOrEmpty(source)){
                    source = "My";
                }
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
        MyCoreBatchRegistrationDialog dialog = new MyCoreBatchRegistrationDialog();
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
}
