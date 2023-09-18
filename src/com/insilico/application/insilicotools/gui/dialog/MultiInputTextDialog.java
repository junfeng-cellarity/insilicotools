package com.insilico.application.insilicotools.gui.dialog;

import com.jidesoft.range.IntegerRange;
import com.jidesoft.range.NumericRange;
import com.jidesoft.range.Range;

import javax.swing.*;
import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.HashMap;

/**
 * Created by jfeng1 on 8/1/16.
 */
public class MultiInputTextDialog extends JDialog {
    JButton okBtn;
    JButton cancelBtn;
    String tag;
    String value;
    boolean isCommitted = false;
    JTextField tagField = new JTextField(5);
    JTextField valueField = new JTextField(5);

    public MultiInputTextDialog(JFrame owner) {
        super(owner);
        setModal(true);
        tagField.setPreferredSize(new Dimension(150,30));
        valueField.setPreferredSize(new Dimension(150,30));
        JPanel contentPanel = new JPanel(new BorderLayout());
        JPanel listPanel = new JPanel(new GridLayout(2,1));

        JPanel p1 = new JPanel(new GridLayout(1,2));
        JPanel p1_1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        p1_1.add(new JLabel("Tag:"));
        p1.add(p1_1);
        p1.add(tagField);
        listPanel.add(p1);

        JPanel p2 = new JPanel(new GridLayout(1,2));
        JPanel p2_1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
        p2_1.add(new JLabel("Value:"));
        p2.add(p2_1);
        p2.add(valueField);
        listPanel.add(p2);


        listPanel.setBorder(BorderFactory.createTitledBorder("Options"));
        contentPanel.add(new JScrollPane(listPanel),BorderLayout.CENTER);
        JPanel btnPanel = new JPanel();
        okBtn = new JButton("OK");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = true;
                setVisible(false);
            }
        });
        btnPanel.add(okBtn);

        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                setVisible(false);
            }
        });
        btnPanel.add(cancelBtn);
        contentPanel.add(btnPanel,BorderLayout.SOUTH);

        setContentPane(contentPanel);
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public String getTag() {
        tag = tagField.getText().trim();
        return tag;
    }

    public String getValue() {
        value = valueField.getText().trim();
        return value;
    }

    public static void main(String[] args) {
        MultiInputTextDialog dialog = new MultiInputTextDialog(null);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setSize(new Dimension(300,180));
        dialog.setVisible(true);

    }
}
