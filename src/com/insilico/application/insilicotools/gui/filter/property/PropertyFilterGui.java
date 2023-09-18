package com.insilico.application.insilicotools.gui.filter.property;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.widget.FloatTextField;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;

public class PropertyFilterGui extends FilterGui {
    private JPanel box;
    private final ArrayList<Line> lines = new ArrayList<Line>();
    private final GridBagConstraints constraints = new GridBagConstraints();
    private final JLabel verticalFillLabel = new JLabel();

    private final PropertyFilterState state;

    public PropertyFilterGui(PropertyFilterState state) {
        this.state = state;

        constraints.gridx = 0;
        constraints.gridy = 99;
        constraints.insets = new Insets(0, 0, 0, 0);
        constraints.weighty = 1.0;
        constraints.fill = GridBagConstraints.VERTICAL;

        final GridBagConstraints labelConstraints = new GridBagConstraints();
        labelConstraints.gridx = 0;
        labelConstraints.gridy = 0;
        labelConstraints.insets = new Insets(0, 0, 0, 0);
        labelConstraints.anchor = GridBagConstraints.NORTHEAST;
        labelConstraints.fill = GridBagConstraints.NONE;
        box = new JPanel(new GridBagLayout());
        box.add(verticalFillLabel, constraints);
        box.setPreferredSize(new Dimension(800, 400));

        PropertyCondition[] conditions = state.getConditions();
        if (conditions != null && conditions.length > 0) {
            lines.clear();
            box.removeAll();
            for (int i = 0; i < conditions.length; i++) {
                Line line = new Line();
                lines.add(line);
                line.textField2.setText(conditions[i].getText2());
                line.textField.setText(conditions[i].getText());
                line.operators.setSelectedItem(conditions[i].getOperator());
                line.properties.setSelectedItem(conditions[i].getProperty());

                if (i == 0) {
                    box.add(buildLinePanel(line, true), labelConstraints);
                    box.add(verticalFillLabel, constraints);
                } else {
                    GridBagConstraints lc = new GridBagConstraints();
                    lc.gridx = 0;
                    lc.gridy = GridBagConstraints.RELATIVE;
                    lc.insets = new Insets(0, 0, 0, 0);
                    lc.anchor = GridBagConstraints.NORTHEAST;
                    lc.fill = GridBagConstraints.NONE;

                    if (lines.size() == 2) lines.get(0).removeButton.setEnabled(true);

                    lines.get(lines.size() - 2).addButton.setEnabled(false);
                    line.addButton.setEnabled(true);

                    box.remove(verticalFillLabel);
                    box.add(buildLinePanel(line, false), lc);
                    box.add(verticalFillLabel, constraints);
                    if (box.getParent() != null) {
                        box.invalidate();
                        box.getParent().validate();
                        box.getParent().repaint();
                    }
                }
            }
        } else {
            Line line = new Line();
            lines.add(line);
            box.add(buildLinePanel(line, true), labelConstraints);
        }

        add(box, BorderLayout.CENTER);
    }

    private JPanel buildLinePanel(final Line l, boolean first) {
        final JPanel p = new JPanel();
        p.add(l.properties);
        p.add(l.operators);
        p.add(l.textField);
        if (l.operators.getSelectedItem().equals(PropertyCondition.BETWEEN)) {
            p.add(l.textField2);
        }
        l.operators.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (l.operators.getSelectedItem().equals(PropertyCondition.BETWEEN)) {
                    p.removeAll();
                    p.add(l.properties);
                    p.add(l.operators);
                    p.add(l.textField);
                    p.add(new JLabel(" and "));
                    p.add(l.textField2);
                    p.add(l.addButton);
                    p.add(l.removeButton);
                } else {
                    p.removeAll();
                    p.add(l.properties);
                    p.add(l.operators);
                    p.add(l.textField);
                    p.add(l.addButton);
                    p.add(l.removeButton);
                }
                if (box.getParent() != null) {
                    box.invalidate();
                    box.getParent().validate();
                    box.getParent().repaint();
                }

            }
        });

        l.addButton.setBorder(null);
        l.addButton.setBorderPainted(false);
        l.addButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Line line = new Line();
                lines.add(line);

                GridBagConstraints labelConstraints = new GridBagConstraints();
                labelConstraints.gridx = 0;
                labelConstraints.gridy = GridBagConstraints.RELATIVE;
                labelConstraints.insets = new Insets(0, 0, 0, 0);
                labelConstraints.anchor = GridBagConstraints.NORTHEAST;
                labelConstraints.fill = GridBagConstraints.NONE;

                if (lines.size() == 2) lines.get(0).removeButton.setEnabled(true);

                lines.get(lines.size() - 2).addButton.setEnabled(false);
                line.addButton.setEnabled(true);

                box.remove(verticalFillLabel);
                box.add(buildLinePanel(line, false), labelConstraints);
                box.add(verticalFillLabel, constraints);
                if (box.getParent() != null) {
                    box.invalidate();
                    box.getParent().validate();
                    box.getParent().repaint();
                }
            }
        });
        l.removeButton.setBorder(null);
        l.removeButton.setBorderPainted(false);
        l.removeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                box.remove(p);
                lines.remove(l);
                if (lines.size() == 1) {
                    lines.get(0).removeButton.setEnabled(false);
                }
                lines.get(lines.size() - 1).addButton.setEnabled(true);

                if (box.getParent() != null) {
                    box.invalidate();
                    box.getParent().validate();
                    box.getParent().repaint();
                }
            }
        });

        if (first) {
            l.removeButton.setEnabled(false);
        }

        p.add(l.addButton);
        p.add(l.removeButton);

        return p;
    }

    public FilterState getState() {
        ArrayList<PropertyCondition> l = new ArrayList<PropertyCondition>();

        for (final Line line : lines) {
            String text = line.textField.getText();
            String property = (String) line.properties.getSelectedItem();
            String operator = (String) line.operators.getSelectedItem();
            String text2 = line.textField2.getText();
            l.add(new PropertyCondition(property, operator, text, text2));
        }
        state.setConditions(l.toArray(new PropertyCondition[l.size()]));
        return state;
    }

    class Line {
        public final FloatTextField textField = new FloatTextField(10);
        public final FloatTextField textField2 = new FloatTextField(10);
        public final JComboBox properties = new JComboBox(PropertyFilter.PROPERTIES_USED_FOR_FILTER);
        public final JComboBox operators = new JComboBox(new String[]{PropertyCondition.EQUALS, PropertyCondition.LESS_THAN, PropertyCondition.GREATER_THAN, PropertyCondition.BETWEEN});
        public final JButton addButton = new JButton(new ImageIcon(getClass().getClassLoader().getResource("add.png")));
        public final JButton removeButton = new JButton(new ImageIcon(getClass().getClassLoader().getResource("remove.png")));
    }
}