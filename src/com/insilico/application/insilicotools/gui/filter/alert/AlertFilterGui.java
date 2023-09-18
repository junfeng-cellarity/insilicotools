package com.insilico.application.insilicotools.gui.filter.alert;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import org.jdesktop.swingx.JXTable;

import javax.swing.*;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.*;
import java.util.List;

public class AlertFilterGui extends FilterGui {
    private final AlertFilterTableModel model;
    private AlertFilterState state;
    private int type = AlertFilter.PAINS_FILTER;

    public AlertFilterGui(AlertFilterState state) {
        this.state = state;

        model = new AlertFilterTableModel(AlertFilterFactory.getAlertRules(AlertFilter.LIBRARY_GUIDELINE_FILTER));
        refresh();

        JXTable ruleTable = new JXTable(model);
        ruleTable.setRowHeight(20);
        ruleTable.getSelectionModel().setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        ruleTable.getColumnModel().getColumn(3).setCellEditor(new SpinnerEditor(0, 0, 100, 1));
        ruleTable.getColumnModel().getColumn(3).setCellRenderer(new SpinnerRenderer(0, 100, 1));

        ruleTable.getColumnModel().getColumn(1).setPreferredWidth(300);

        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        add(new JLabel("Select the structure alert rules you want to filter on."));
        add(new JScrollPane(ruleTable));

        JToolBar toolbar = new JToolBar("Compounds", JToolBar.HORIZONTAL);
        toolbar.setFloatable(false);

        final JComboBox cb = new JComboBox(new String[]{"Library","PAINS"});
        cb.setMaximumSize(new Dimension(150,20));
        cb.setMinimumSize(new Dimension(150,20));

        cb.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String ruleName = (String)cb.getSelectedItem();
                AlertFilterGui.this.firePropertyChange("nameChanged","",ruleName);
                if(ruleName.equals("PAINS")){
                    type = AlertFilter.PAINS_FILTER;
                }else{
                    type = AlertFilter.LIBRARY_GUIDELINE_FILTER;
                }
                model.setRules(AlertFilterFactory.getAlertRules(type));
            }
        });
        toolbar.add(new JLabel("Rule Set:"));
        toolbar.add(cb);
        toolbar.addSeparator();
        JButton selectAllButton = new JButton("Select All");
        selectAllButton.setBorder(null);
        selectAllButton.setBorderPainted(false);
        selectAllButton.setToolTipText("Select All");
        selectAllButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                for (int i = 0; i < model.getRowCount(); i++) {
                    model.setValueAt(Boolean.TRUE, i, 2);
                }
                model.fireTableDataChanged();
            }
        });
        toolbar.add("Select All", selectAllButton);
        toolbar.addSeparator();

        JButton clearSelectionButton = new JButton("Clear Selection");
        clearSelectionButton.setBorder(null);
        clearSelectionButton.setBorderPainted(false);
        clearSelectionButton.setToolTipText("Clear Selection");
        clearSelectionButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                for (int i = 0; i < model.getRowCount(); i++) {
                    model.setValueAt(Boolean.FALSE, i, 2);
                }
                model.fireTableDataChanged();
            }
        });
        toolbar.add("Clear Selection", clearSelectionButton);

        add(toolbar);
    }

    public void refresh() {
        Vector<AlertRule> selectedRules = state.getSelectedRules();

        for (int i = 0; i < model.getRowCount(); i++) {
            if (selectedRules.size() == 0) {
                model.setValueAt(Boolean.TRUE, i, 2);
            } else {
                model.setValueAt(Boolean.FALSE, i, 2);
                for (int j = 0; j < selectedRules.size(); j++) {
                    if (selectedRules.get(j).getIndex() == model.getRuleAt(i).getIndex()) {
                        model.setValueAt(Boolean.TRUE, i, 2);
                        model.setValueAt(new Integer(selectedRules.get(j).getMax()), i, 3);
                        break;
                    }
                }
            }
        }
        model.fireTableDataChanged();
    }

    public FilterState getState() {
        Vector<AlertRule> selectedRules = new Vector<AlertRule>();
        for (int i = 0; i < model.getRowCount(); i++) {
            if (((Boolean) model.getValueAt(i, 2)).booleanValue()) selectedRules.add((AlertRule) model.getValueAt(i, 1));
        }
        state.setSelectedRules(selectedRules);

        return state;
    }

    class AlertFilterTableModel extends AbstractTableModel{
        private List<AlertRule> rules;

        public AlertFilterTableModel(List<AlertRule> rules) {
            this.rules = rules;
        }

        public void setRules(List<AlertRule> rules) {
            this.rules = rules;
            fireTableStructureChanged();
        }

        public int getColumnCount() {
            return 4;
        }

        public int getRowCount() {
            return rules.size();
        }

        public AlertRule getRuleAt(int rowIndex) {
            return (AlertRule) rules.get(rowIndex);
        }

        public Object getValueAt(int rowIndex, int columnIndex) {
            AlertRule rule = (AlertRule) rules.get(rowIndex);
            switch (columnIndex) {
                case 0:
                    return new Integer(rule.getIndex());

                case 1:
                    return rule;

                case 2:
                    return Boolean.valueOf(rule.isSelected());

                case 3:
                    return new Integer(rule.getMax());

                default:
                    throw new IllegalArgumentException("Unknown switch state");
            }
        }

        public boolean isCellEditable(int rowIndex, int columnIndex) {
            return columnIndex > 1;
        }

        public Class<?> getColumnClass(int columnIndex) {
            switch (columnIndex) {
                case 0:
                    return Integer.class;

                case 1:
                    return String.class;

                case 2:
                    return Boolean.class;

                case 3:
                    return Integer.class;

                default:
                    throw new IllegalArgumentException("Unknown switch state");
            }
        }

        public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
            AlertRule rule = (AlertRule) rules.get(rowIndex);
            if (columnIndex == 3) {
                rule.setMax(Integer.parseInt(aValue.toString()));
            } else if (columnIndex == 2) {
                rule.setSelected(((Boolean) aValue).booleanValue());
            }
        }

        public String getColumnName(int column) {
            switch (column) {
                case 0:
                    return "Rule ID";

                case 1:
                    return "Reos Rule";

                case 2:
                    return "Selected";

                case 3:
                    return "Count";

                default:
                    throw new IllegalArgumentException("Unknown switch state");
            }
        }

        public void reset() {
            rules = AlertFilterFactory.getAlertRules(type);
        }

        public void update(int id, int count) {
            Iterator<AlertRule> itr = rules.iterator();
            while (itr.hasNext()) {
                AlertRule holder = (AlertRule) itr.next();
                if (holder.getIndex() == id) {
                    holder.setSelected(true);
                    holder.setMax(count);
                    break;
                }
            }
            fireTableDataChanged();
        }

        public boolean isSortable(int column) {
            return true;
        }

        public void sortColumn(final int column, final boolean ascending) {
            Collections.sort(rules, new Comparator<AlertRule>() {
                public int compare(AlertRule o1, AlertRule o2) {
                    int result = 0;

                    switch (column) {
                        case 0:
                            result = new Integer(((AlertRule) o1).getIndex()).compareTo(new Integer(((AlertRule) o2).getIndex()));
                            break;

                        case 1:
                            result = ((AlertRule) o1).getDescription().compareToIgnoreCase(((AlertRule) o2).getDescription());
                            break;

                        case 2:
                            if (((AlertRule) o1).isSelected() == ((AlertRule) o2).isSelected()) {
                                result = 0;
                            } else if (((AlertRule) o1).isSelected()) {
                                result = 1;
                            } else {
                                result = -1;
                            }
                            break;

                        case 3:
                            result = new Integer(o1.getMax()).compareTo(new Integer(o2.getMax()));
                            break;
                    }

                    return ascending ? result : 1 - result;
                }
            });
        }

        public AlertRule[] getRules() {
            return rules.toArray(new AlertRule[rules.size()]);
        }
    }

    class SpinnerEditor extends AbstractCellEditor implements TableCellEditor {
        final JSpinner spinner = new JSpinner();

        public SpinnerEditor(int value, int min, int max, int step) {
            spinner.setModel(new SpinnerNumberModel(value, min, max, step));
            spinner.setBorder(null);
        }

        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
            spinner.setValue(value);
            model.setValueAt(Boolean.TRUE, row, 1);
            return spinner;
        }

        public boolean isCellEditable(EventObject evt) {
            return true;
        }

        public Object getCellEditorValue() {
            return spinner.getValue();
        }
    }

    class SpinnerRenderer extends JSpinner implements TableCellRenderer {
        public SpinnerRenderer(int min, int max, int step) {
            super(new SpinnerNumberModel(0, min, max, step));
            setBorder(null);
        }

        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
            setValue(value);

            if (isSelected) {
                setBackground(table.getSelectionBackground());
            } else {
                setBackground(table.getBackground());
            }
            return this;
        }
    }
}