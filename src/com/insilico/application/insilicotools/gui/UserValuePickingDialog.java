package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.util.ChemFunc;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.awt.event.*;
import java.util.Arrays;
import java.util.Vector;

/**
 * Created by jfeng1 on 10/13/15.
 */
public class UserValuePickingDialog extends JDialog {
    boolean isSubmitted = false;
    JComboBox tagCB;

    String selectedTag;


    public UserValuePickingDialog(Vector<String> userTags) {
        JPanel contentPane = new JPanel(new BorderLayout());
        JPanel gridPanel;
        tagCB = new JComboBox(userTags);
        gridPanel = new JPanel();
        gridPanel.add(new JLabel("Select NameTag:"));
        gridPanel.add(tagCB);
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
            setSize(new Dimension(400,120));
    }

    private void onOK(){
        isSubmitted = true;
        selectedTag = (String)tagCB.getSelectedItem();
        setVisible(false);
    }

    private void onCancel(){
        isSubmitted = false;
        setVisible(false);
    }

    public String getSelectedTag() {
        return selectedTag;
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
        UserValuePickingDialog dialog = new UserValuePickingDialog(new Vector<String>(Arrays.asList(properties)));
        dialog.setVisible(true);
    }


}
