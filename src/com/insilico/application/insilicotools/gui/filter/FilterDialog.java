package com.insilico.application.insilicotools.gui.filter;

import com.insilico.application.insilicotools.data.PropertyMolecule;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Created by jfeng1 on 4/13/16.
 */
public class FilterDialog extends JDialog{
    TreeFilter filterPanel;
    PropertyMolecule[] selectedMolecules;
    boolean isCommitted = false;

    public FilterDialog(JFrame owner, PropertyMolecule[] propertyMolecules) {
        super(owner, "Compound Filtering");
        initialize(propertyMolecules);
    }

    private void initialize(PropertyMolecule[] propertyMolecules) {
        setModal(true);
        setGlassPane(new GhostGlassPane());
        JPanel p = new JPanel(new BorderLayout());
        final JButton okBtn = new JButton("OK");
        filterPanel = new TreeFilter(propertyMolecules);
        filterPanel.addShoppingCartListener(new MoleculeSelectionListener() {
            @Override
            public void moleculesSelected(MoleculeSelectionEvent evt) {
                selectedMolecules = evt.getMols();
                if(selectedMolecules.length>0){
                    okBtn.setEnabled(true);
                }else{
                    okBtn.setEnabled(false);
                }
            }
        });
        JPanel btnPanel = new JPanel();
        okBtn.setEnabled(false);
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = true;
                FilterDialog.this.setVisible(false);
            }
        });


        JButton cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                FilterDialog.this.setVisible(false);
            }
        });
        btnPanel.add(okBtn);
        btnPanel.add(cancelBtn);
        p.add(filterPanel,BorderLayout.CENTER);
        p.add(btnPanel,BorderLayout.SOUTH);
        getContentPane().add(p);
        setSize(new Dimension(800,600));
    }

    public FilterDialog(JDialog owner, PropertyMolecule[] propertyMolecules) {
        super(owner, "Compound Filtering");
        initialize(propertyMolecules);
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public PropertyMolecule[] getSelectedMolecules() {
        return selectedMolecules;
    }

}
