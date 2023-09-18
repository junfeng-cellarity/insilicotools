package com.insilico.application.insilicotools.gui.filter;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;

import javax.swing.*;
import java.awt.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

public abstract class GuiFilterController extends FilterController {
    private FilterGui gui;

    protected GuiFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected abstract FilterGui getGuiComponent(FilterState state);

    public int showOptionDialog(JComponent parentComponent, String title, Object message, int messageType, int optionType, Icon icon, Object[] options, Object initialValue) {
        JOptionPane pane = new JOptionPane(message, messageType, optionType, icon, options, initialValue);

        JDialog dialog = pane.createDialog(parentComponent, title);
        dialog.setResizable(true);
        dialog.setVisible(true);
        dialog.pack();
        dialog.dispose();

        Object selectedValue = pane.getValue();

        if(selectedValue == null) return JOptionPane.CLOSED_OPTION;
        if(options == null) {
            if(selectedValue instanceof Integer) return ((Integer) selectedValue).intValue();
        }
        else {
            for(int counter = 0, maxCounter = options.length; counter < maxCounter; counter++) {
                if(options[counter].equals(selectedValue)) return counter;
            }
        }
        return JOptionPane.CLOSED_OPTION;
    }

    public FilterResult filter(ProgressReporter reporter, PropertyMolecule[] mols) throws Exception {
        if (state.isNewFilter()) {
            if (gui == null) gui = getGuiComponent(state);

            JPanel p = new JPanel(new BorderLayout());
            p.add(gui, BorderLayout.CENTER);
            JPanel namePanel = new JPanel();
            p.add(namePanel, BorderLayout.SOUTH);
            final JTextField nameField = new JTextField(getName(), 20);
            gui.addPropertyChangeListener(new PropertyChangeListener() {
                @Override
                public void propertyChange(PropertyChangeEvent evt) {
                    if(evt.getPropertyName().equals("nameChanged")){
                        nameField.setText((String)evt.getNewValue());
                    }
                }
            });
            namePanel.add(new JLabel("Name:"));
            namePanel.add(nameField);

            while (true) {
                gui.onScreen();
                String name = showOptionDialog(treeFilter, getName(), p, JOptionPane.PLAIN_MESSAGE, JOptionPane.OK_CANCEL_OPTION, null, null, null) == JOptionPane.OK_OPTION ? nameField.getText() : null;
                gui.offScreen();
                if (name == null) return null;
                String valid = gui.isValidState();
                if (valid == null) {
                    setName(name);
                    state = gui.getState();
                    state.setNewFilter(false);
                    return super.filter(reporter, mols);
                } else {
                    JOptionPane.showMessageDialog(null, valid, "Warning", JOptionPane.PLAIN_MESSAGE);
                }
            }
        } else{
            return super.filter(reporter, mols);
        }
    }
}