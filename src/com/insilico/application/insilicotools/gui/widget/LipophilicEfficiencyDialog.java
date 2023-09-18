package com.insilico.application.insilicotools.gui.widget;

import com.insilico.application.insilicotools.util.ChemFunc;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.Arrays;
import java.util.Vector;

/**
 * Created by jfeng1 on 6/6/16.
 */
public class LipophilicEfficiencyDialog extends JDialog {
    boolean isSubmitted = false;
    JComboBox tagCB;
    JComboBox logpTypeCB;
    JComboBox unitCB;

    String selectedTag;
    String unit;
    String logpType;
//    boolean useLog;


    public LipophilicEfficiencyDialog(Vector<String> userTags) {
        JPanel contentPane = new JPanel(new BorderLayout());
        JPanel gridPanel;
        tagCB = new JComboBox(userTags);
        gridPanel = new JPanel(new GridLayout(3,1));
        unitCB = new JComboBox(new String[]{"uM", "nM"});
        logpTypeCB = new JComboBox(new String[]{"CLogP","MoKa LogP","XLogP","ChemAxon LogP"});
//        JCheckBox useLogCB = new JCheckBox("-Log",false);

        JPanel p1 = new JPanel();
        p1.add(new JLabel("Assay Name:"));
        p1.add(tagCB);

        JPanel p2 = new JPanel();
        p2.add(new JLabel("Unit:"));
        p2.add(unitCB);

        JPanel p3 = new JPanel();
        p3.add(new JLabel("LogP Type:"));
        p3.add(logpTypeCB);
//        p2.add(useLogCB);

        gridPanel.add(p1);
        gridPanel.add(p2);
        gridPanel.add(p3);
        gridPanel.setBorder(new EtchedBorder());

        contentPane.add(gridPanel, BorderLayout.CENTER);

        JButton buttonOK = new JButton("OK");
        buttonOK.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onOK();
            }
        });

        JButton buttonCancel = new JButton("Cancel");
        buttonCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        });

        JPanel btnPanel = new JPanel();
        btnPanel.add(buttonOK);
        btnPanel.add(buttonCancel);
        contentPane.add(btnPanel, BorderLayout.PAGE_END);
        setContentPane(contentPane);

// call onCancel() when cross is clicked
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setModal(true);
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                onCancel();
            }
        });

// call onCancel() on ESCAPE
        contentPane.registerKeyboardAction(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        }, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
        setSize(new Dimension(400, 240));
    }

    private void onOK(){
        isSubmitted = true;
        selectedTag = (String)tagCB.getSelectedItem();
        unit = (String) unitCB.getSelectedItem();
        logpType = (String)logpTypeCB.getSelectedItem();
        setVisible(false);
    }

    private void onCancel(){
        isSubmitted = false;
        setVisible(false);
    }

    public String getSelectedTag() {
        return selectedTag;
    }

    public String getUnit() {
        return unit;
    }

    public String getLogpType() {
        return logpType;
    }

    public boolean isSubmitted() {
        return isSubmitted;
    }

    @Override
    public void setVisible(boolean b) {
        if(b){
            isSubmitted=false;
        }
        super.setVisible(b);
    }

    public static void main(String[] args) {
        String[] properties = ChemFunc.properties;
        LipophilicEfficiencyDialog dialog = new LipophilicEfficiencyDialog(new Vector<String>(Arrays.asList(properties)));
        dialog.setVisible(true);
    }


}
