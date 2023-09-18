package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.PropertyMolecule;

import javax.swing.table.DefaultTableModel;
import java.util.List;
import java.util.Vector;

/**
 * Created by jfeng1 on 12/21/15.
 */
public class MatrixMolTableModel extends DefaultTableModel{
    Vector<PropertyMolecule> propertyMolecules;
    int columnCount = 4;

    public Vector<PropertyMolecule> getPropertyMolecules() {
        return propertyMolecules;
    }

    public MatrixMolTableModel(Vector<PropertyMolecule> molecules) {
        propertyMolecules = molecules;
        if(molecules==null){
            propertyMolecules = new Vector<PropertyMolecule>();
        }
    }

    public void addMolecules(List<PropertyMolecule> mols){
        if(mols==null||mols.isEmpty()){
            return;
        }
        for(PropertyMolecule m:mols){
            if(!propertyMolecules.contains(m)) {
                propertyMolecules.add(m);
            }
        }
        fireTableStructureChanged();
    }

    @Override
    public void setColumnCount(int columnCount) {
        this.columnCount = columnCount;
        this.fireTableStructureChanged();
        this.fireTableDataChanged();
    }

    public void clear(){
        propertyMolecules.clear();
        this.fireTableStructureChanged();
    }

    @Override
    public int getRowCount() {
        if(propertyMolecules==null||propertyMolecules.isEmpty()){
            return 0;
        }
        int a = (propertyMolecules.size()%columnCount)==0?0:1;
        return propertyMolecules.size()/columnCount+a;
    }

    @Override
    public int getColumnCount() {
        return columnCount;
    }

    @Override
    public String getColumnName(int column) {
        return ""+column%columnCount;
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        return false;
    }

    @Override
    public void setValueAt(Object aValue, int row, int column) {
    }

    @Override
    public Object getValueAt(int row, int column) {
        if(propertyMolecules.size()==0){
            return null;
        }
        if(propertyMolecules.size()<=row*columnCount+column){
            return null;
        }
        return  propertyMolecules.get(row * columnCount + column);
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        return PropertyMolecule.class;
    }

}
