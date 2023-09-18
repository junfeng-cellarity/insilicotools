package com.insilico.application.insilicotools.gui.dialog;

import com.insilico.application.insilicotools.util.ChemFunc;
import com.jidesoft.swing.CheckBoxList;
import com.jidesoft.swing.JideScrollPane;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

/**
 * Created by jfeng1 on 3/17/16.
 */
public class MultiChoicesDialog extends JDialog {
    boolean isSubmitted = false;
    Vector<String> selectedTags;
    CheckBoxList checkBoxList;
//    boolean useLog;

    public MultiChoicesDialog(JFrame frame,String[] userTags, Vector<String> existingTags) {
        super(frame);
        final JPanel contentPane = new JPanel(new BorderLayout());
        checkBoxList = new CheckBoxList(userTags);
        selectedTags = new Vector<String>();
        contentPane.add(new JideScrollPane(checkBoxList));
        if(existingTags!=null&&!existingTags.isEmpty()) {
            selectedTags.addAll(existingTags);
            Vector<Integer> indices = new Vector<Integer>();
            for (int i=0;i<userTags.length;i++) {
                String tag = userTags[i];
                if(existingTags.contains(tag)){
                    indices.add(i);
                }
            }
            int[] indices2 = new int[indices.size()];
            for(int i=0;i<indices.size();i++){
                indices2[i] = indices.get(i);
            }
            checkBoxList.setCheckBoxListSelectedIndices(indices2);
        }

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
        setSize(new Dimension(400,400));
        setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setModal(true);
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                onCancel();
            }
        });

        final JPopupMenu menu = new JPopupMenu();
        JMenuItem selectallItem = new JMenuItem("Select all");
        selectallItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                checkBoxList.selectAll();
            }
        });
        menu.add(selectallItem);

        JMenuItem unselectAllItem = new JMenuItem("Select none");
        unselectAllItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                checkBoxList.selectNone();
            }
        });
        menu.add(unselectAllItem);

        checkBoxList.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if(SwingUtilities.isRightMouseButton(e)) {
                    menu.show(checkBoxList, e.getX(), e.getY());
                }
            }
        });

// call onCancel() on ESCAPE
        contentPane.registerKeyboardAction(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                onCancel();
            }
        }, KeyStroke.getKeyStroke(KeyEvent.VK_ESCAPE, 0), JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT);
    }

    public boolean isSubmitted() {
        return isSubmitted;
    }

    private void onOK(){
        isSubmitted = true;
        selectedTags.clear();
        Object[] selectedList = checkBoxList.getCheckBoxListSelectedValues();
        for(Object value:selectedList){
            selectedTags.add((String)value);
        }
        setVisible(false);
    }

    private void onCancel(){
        isSubmitted = false;
        setVisible(false);
    }

    public Vector<String> getSelectedTags() {
        return selectedTags;
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
        MultiChoicesDialog dialog = new MultiChoicesDialog(null,properties,null);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setVisible(true);
        System.out.println(dialog.getSelectedTags());
    }


}
