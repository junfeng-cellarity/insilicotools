package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.util.ChemFunc;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

public class PropertySelectionDialog extends JDialog {
    private PropertySelectionPanel contentPane;
    private JButton buttonOK;
    private JButton buttonCancel;
    boolean isSubmitted = false;

    @Override
    public PropertySelectionPanel getContentPane() {
        return contentPane;
    }


    public PropertySelectionDialog(String[] properties, String[] adme_models, String[] safety_models, String[] more_models) {
        contentPane = new PropertySelectionPanel(properties, adme_models,safety_models,more_models);
        setContentPane(contentPane);
        setModal(true);
        getRootPane().setDefaultButton(buttonOK);
        buttonOK = new JButton("OK");
        buttonOK.setEnabled(true);
        buttonCancel = new JButton("Cancel");

        buttonOK.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onOK();
            }
        });

        buttonCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        });
        JPanel btnPanel = new JPanel();
        btnPanel.add(buttonOK);
        btnPanel.add(buttonCancel);
        contentPane.addPropertyChangeListener("selection", new PropertyChangeListener() {
            @Override
            public void propertyChange(PropertyChangeEvent evt) {
                if(contentPane.getSelectedProperties().isEmpty()){
                    buttonOK.setEnabled(false);
                }else{
                    buttonOK.setEnabled(true);
                }
            }
        });
        contentPane.add(btnPanel, BorderLayout.PAGE_END);

// call onCancel() when cross is clicked
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
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
        setSize(new Dimension(800, 515));
    }

    private void onOK() {
        isSubmitted = true;
        setVisible(false);
    }

    private void onCancel() {
// add your code here if necessary
        isSubmitted = false;
        setVisible(false);
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
        String[] models = ChemFunc.adme_models;
        String[] safety_models = ChemFunc.safety_models;
        String[] other_modesl = ChemFunc.other_models;
        PropertySelectionDialog dialog = new PropertySelectionDialog(properties,models, safety_models,other_modesl);
        dialog.setSize(800, 600);
        dialog.setVisible(true);
        System.exit(0);
    }
}
