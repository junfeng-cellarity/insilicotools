package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.SVGTableCellRenderer;
import org.jdesktop.swingx.JXTable;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.DefaultTableModel;
import java.awt.*;

/**
 * Created by jfeng1 on 5/13/16.
 */
public class MolTable2D extends JXTable{
    DefaultTableModel tableModel;
    int molTableRowHeight = 150;

    public MolTable2D(DefaultTableModel tableModel) {
        super(tableModel);
        this.tableModel = tableModel;
        setShowGrid(true);
        setGridColor(Color.GRAY);
        setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer());
        setRowHeight(300);
        setSortOrderCycle(SortOrder.ASCENDING, SortOrder.DESCENDING, SortOrder.UNSORTED);
        setColumnControlVisible(true);
        this.getColumnModel().addColumnModelListener(new TableColumnModelListener() {
            @Override
            public void columnAdded(TableColumnModelEvent e) {

            }

            @Override
            public void columnRemoved(TableColumnModelEvent e) {

            }

            @Override
            public void columnMoved(TableColumnModelEvent e) {

            }

            @Override
            public void columnMarginChanged(ChangeEvent e) {
                if (MolTable2D.this.getColumnModel().getColumnCount() > 0) {
                    int molColIdx = MolTable2D.this.convertColumnIndexToModel(1);
                    molTableRowHeight = MolTable2D.this.getColumnModel().getColumn(molColIdx).getWidth();
                    MolTable2D.this.setRowHeight(molTableRowHeight);
                }
            }

            @Override
            public void columnSelectionChanged(ListSelectionEvent e) {

            }
        });

    }

    public void updateTable(){
        tableModel.fireTableDataChanged();
        tableModel.fireTableStructureChanged();
        fixTableFormat();
    }

    public void fixTableFormat(){
        if (tableModel.getColumnCount() > 0) {
            setRowHeight(molTableRowHeight);
            setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
            int selectColIdx = convertColumnIndexToModel(0);
            int molColIdx = convertColumnIndexToModel(1);
            getColumnModel().getColumn(selectColIdx).setMaxWidth(20);
            getColumnModel().getColumn(molColIdx).setMinWidth(150);
        }
    }

}
