package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.SVGTableCellRenderer;
import org.jdesktop.swingx.JXTable;

import javax.swing.*;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.TableModelEvent;
import javax.swing.event.TableModelListener;
import java.util.ArrayList;
import java.util.Vector;

/**
 * Created by jfeng1 on 2/17/16.
 */

class Cell {
    private int row;

    private int column;

    public Cell(int row, int column) {
        this.row = row;
        this.column = column;
    }

    public boolean is(int r, int c) {
        return row == r && column == c;
    }
}

class CellSelectionSet {
    private ArrayList<Cell> cells = new ArrayList<Cell>();

    public int size(){
        return cells.size();
    }

    public void add(int r, int c) {
        if (!contains(r, c)) {
            cells.add(new Cell(r, c));
        }
    }

    public void remove(int r, int c){
        Cell cellToRemove = null;
        if(contains(r,c)){
            for(Cell cell:cells){
                if(cell.is(r,c)){
                    cellToRemove = cell;
                    break;
                }
            }
        }
        if(cellToRemove!=null){
            cells.remove(cellToRemove);
        }
    }

//    public boolean containsOneOrLess() {
//        return cells.size() <= 1;
//    }

    public boolean contains(int r, int c) {
        for (Cell cell : cells) {
            if (cell.is(r, c)) {
                return true;
            }
        }
        return false;
    }

    public void clear() {
        cells.clear();
    }


}

public class MatrixCellTable extends JXTable{
    CellSelectionSet cellSelectionSet = new CellSelectionSet();

    @Override
    public boolean isCellEditable(int row, int column) {
        return false;
    }

    @Override
    public void changeSelection(int rowIndex, int columnIndex, boolean toggle, boolean extend) {
        super.changeSelection(rowIndex, columnIndex, toggle, extend);
        if(cellSelectionSet.contains(rowIndex,columnIndex)){
            if(!extend){
                cellSelectionSet.remove(rowIndex,columnIndex);
            }
            return;
        }
        if (toggle) {
            cellSelectionSet.add(rowIndex, columnIndex);

        } else {
            if (extend) {
                cellSelectionSet.add(rowIndex, columnIndex);

            } else {
                // reset
                cellSelectionSet.clear();
                cellSelectionSet.add(rowIndex, columnIndex);
            }
        }

    }

    public void selectAll(){
        MatrixMolTableModel dm = (MatrixMolTableModel)getModel();
        for(int i=0;i<dm.getRowCount();i++){
            for(int j=0;j<dm.getColumnCount();j++){
                Object value = dm.getValueAt(i, j);
                if(value!=null&&value instanceof PropertyMolecule){
                    if(!cellSelectionSet.contains(i,j)){
                        cellSelectionSet.add(i,j);
                        ((PropertyMolecule)value).setIsSelected(true);
                    }
                }
            }
        }
        dm.fireTableDataChanged();
    }

    public void unselectAll(){
        MatrixMolTableModel dm = (MatrixMolTableModel)getModel();
        for(int i=0;i<dm.getRowCount();i++){
            for(int j=0;j<dm.getColumnCount();j++){
                Object value = dm.getValueAt(i, j);
                if(value!=null&&value instanceof PropertyMolecule){
                    if(cellSelectionSet.contains(i,j)){
                        cellSelectionSet.remove(i,j);
                        ((PropertyMolecule)value).setIsSelected(false);
                    }
                }
            }
        }
        dm.fireTableDataChanged();
    }

    public void invertSelection(){
        MatrixMolTableModel dm = (MatrixMolTableModel)getModel();
        for(int i=0;i<dm.getRowCount();i++){
            for(int j=0;j<dm.getColumnCount();j++){
                Object value = dm.getValueAt(i, j);
                if(value!=null&&value instanceof PropertyMolecule) {
                    if (!cellSelectionSet.contains(i, j)) {
                        cellSelectionSet.add(i, j);
                        ((PropertyMolecule)value).setIsSelected(true);
                    } else {
                        cellSelectionSet.remove(i, j);
                        ((PropertyMolecule)value).setIsSelected(false);
                    }
                }
            }
        }
        dm.fireTableDataChanged();
    }

    public int getNumSelected(){
        return cellSelectionSet.size();
    }

    public Vector<PropertyMolecule> getSelectedMols(){
        Vector<PropertyMolecule> selectedMols = new Vector<PropertyMolecule>();
        MatrixMolTableModel dm = (MatrixMolTableModel)getModel();
        for(int i=0;i<dm.getRowCount();i++){
            for(int j=0;j<dm.getColumnCount();j++){
                Object value = dm.getValueAt(i, j);
                if(value!=null&&value instanceof PropertyMolecule) {
                    if (cellSelectionSet.contains(i, j)) {
                        PropertyMolecule selectedMol = (PropertyMolecule) value;
                        selectedMols.add(selectedMol);
                    }
                }
            }
        }
        return selectedMols;
    }

    @Override
    public boolean isCellSelected(int row, int column) {
//        if (cellSelectionSet.containsOneOrLess()) {
//            // show the default
//            return super.isCellSelected(row, column);
//        }
        return cellSelectionSet.contains(row, column);
    }

    public MatrixCellTable(final MatrixMolTableModel dm) {
        super(dm);
        setShowGrid(true);
        setCellSelectionEnabled(true);
        setRowSelectionAllowed(true);
        setColumnSelectionAllowed(true);
        setSortable(false);
        DefaultListSelectionModel selectionModel = new DefaultListSelectionModel();
        selectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        setRowHeight(150);
        setDefaultRenderer(PropertyMolecule.class, new SVGTableCellRenderer(4));
        getSelectionModel().addListSelectionListener(new ListSelectionListener() {
            @Override
            public void valueChanged(ListSelectionEvent e) {
                if(!e.getValueIsAdjusting()){
                    MatrixCellTable table = MatrixCellTable.this;
                    for(int i=0;i<table.getRowCount();i++){
                        for(int j=0;j<table.getColumnCount();j++){
                            Object value = dm.getValueAt(i, j);
                            if(table.isCellSelected(i,j)){
                                if(value!=null&&value instanceof PropertyMolecule){
                                    ((PropertyMolecule)value).setIsSelected(true);
                                }
                            }else{
                                if(value!=null&&value instanceof PropertyMolecule){
                                    ((PropertyMolecule)value).setIsSelected(false);
                                }
                            }
                        }
                    }
                    dm.fireTableDataChanged();
                }
            }
        });

        dm.addTableModelListener(new TableModelListener() {
            @Override
            public void tableChanged(TableModelEvent e) {
                for(int i=0;i<dm.getRowCount();i++){
                    for(int j=0;j<dm.getColumnCount();j++){
                        Object value = dm.getValueAt(i, j);
                        if(value!=null&&value instanceof PropertyMolecule) {
                            PropertyMolecule m = (PropertyMolecule)value;
                            if(m.isSelected()){
                                if (!cellSelectionSet.contains(i, j)) {
                                    cellSelectionSet.add(i, j);
                                }
                            }else{
                                if (cellSelectionSet.contains(i, j)) {
                                    cellSelectionSet.remove(i, j);
                                }
                            }
                        }
                    }
                }
            }
        });



    }
}
