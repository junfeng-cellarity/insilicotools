package com.insilico.application.insilicotools.gui.widget;

import javax.swing.*;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableModel;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.util.EventObject;

class StarRaterCellEditor extends AbstractCellEditor implements TableCellEditor {
    int rank = 0;
    int row = -1;
    StarRater rater = new StarRater(5,rank);
    TableModel tableModel;

    public StarRaterCellEditor() {
        rater.addStarListener(new StarRater.StarListener() {
            @Override
            public void handleSelection(int selection) {
                if(rank!=selection) {
                    rank = selection;
                    if(tableModel!=null){
                        fireEditingStopped();
                        ((DefaultTableModel)tableModel).fireTableRowsUpdated(row,row);
                    }
                }
            }
        });
    }

    @Override
    public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column) {
        if (value instanceof Integer) {
            this.rank = (Integer) value;
        }
        tableModel = table.getModel();
        this.row = row;
        rater.setRating(rank);
        return rater;
    }

    @Override
    public boolean isCellEditable(EventObject e) {
        if (e instanceof MouseEvent) {
            return ((MouseEvent)e).getClickCount() >= 2;
        }
        return true;
    }

    @Override
    public Object getCellEditorValue() {
        return this.rank;
    }
}
