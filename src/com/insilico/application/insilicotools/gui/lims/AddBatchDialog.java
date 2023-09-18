package com.insilico.application.insilicotools.gui.lims;

import com.insilico.application.insilicotools.database.LimsDAO;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.google.common.base.Strings;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.SQLException;

public class AddBatchDialog extends StandardDialog{
    JPanel mainPanel = new JPanel();
    JButton okBtn = new JButton("OK");
    JButton cancelBtn = new JButton("Cancel");
    boolean isCommitted = false;
    JTextField nameField = new JTextField(10);
    Batch batch;
    String userName;

    public AddBatchDialog() {
        super();
        setModal(true);
        initializeComponents();
        setPreferredSize(new Dimension(400,200));
        pack();
    }


    public boolean isCommitted() {
        return isCommitted;
    }

    private void initializeComponents(){
        BoxLayout layout = new BoxLayout(mainPanel,BoxLayout.PAGE_AXIS);
        mainPanel.setLayout(layout);
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Batch Name:"));
        p1.add(nameField);
        JPanel p2 = new JPanel();
        p2.add(new JLabel("Scientist:"));

        try {
            userName = InSlilicoPanel.getInstance().getUserName();
        } catch (Exception e) {
            userName = "Unknown Chemist";
        }
        JTextField dummyField = new JTextField(10);
        dummyField.setEditable(false);
        dummyField.setText(userName);
        p2.add(dummyField);
        mainPanel.add(p1);
        mainPanel.add(p2);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                if(Strings.isNullOrEmpty(nameField.getText())){
                    JOptionPane.showMessageDialog(AddBatchDialog.this,"Name field is empty.");
                    return;
                }
                try {
                    batch = LimsDAO.getInstance().getBatch(LimsDAO.getInstance().insertBatch(nameField.getText().trim(),userName));
                } catch (SQLException e1) {
                    e1.printStackTrace();
                    JOptionPane.showMessageDialog(AddBatchDialog.this,e1.getMessage());
                    return;
                }
                isCommitted = true;
                AddBatchDialog.this.dispose();
            }
        });
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                AddBatchDialog.this.dispose();
            }
        });
    }


    public Batch getBatch() {
        return batch;
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

    public static void main(String[] args) {
        AddBatchDialog dialog = new AddBatchDialog();
        dialog.setVisible(true);
    }
}
