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
public class LigandEfficiencyDialog extends JDialog{
    boolean isSubmitted = false;
    JComboBox tagCB;
    JComboBox unitCB;

    String selectedTag;
    String unit;
//    boolean useLog;


    public LigandEfficiencyDialog(Vector<String> userTags) {
        JPanel contentPane = new JPanel(new BorderLayout());
        JPanel gridPanel;
        tagCB = new JComboBox(userTags);
        gridPanel = new JPanel(new GridLayout(3,1));
        unitCB = new JComboBox(new String[]{"uM", "nM"});

        JPanel p1 = new JPanel();
        p1.add(new JLabel("Assay Name:"));
        p1.add(tagCB);

        JPanel p2 = new JPanel();
        p2.add(new JLabel("Unit:"));
        p2.add(unitCB);


        gridPanel.add(p1);
        gridPanel.add(p2);
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
        setSize(new Dimension(400, 180));
    }

    private void onOK(){
        isSubmitted = true;
        selectedTag = (String)tagCB.getSelectedItem();
        if(unitCB!=null) {
            unit = (String) unitCB.getSelectedItem();
        }
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
        LigandEfficiencyDialog dialog = new LigandEfficiencyDialog(new Vector<String>(Arrays.asList(properties)));
        dialog.setVisible(true);
    }


}
