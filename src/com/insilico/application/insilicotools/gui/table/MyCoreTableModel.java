package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.Core;
import com.insilico.application.insilicotools.data.PropertyMolecule;

import javax.swing.table.DefaultTableModel;
import java.util.Date;
import java.util.Vector;

/**
 * Created by jfeng1 on 2/13/17.
 */
public class MyCoreTableModel extends DefaultTableModel{
    Vector<Core> MyCores;
    String[] columnNames = new String[]{"","Core Structure","Name","Description", "Chemist","Date","Source", "Comment", "Status"};

    public MyCoreTableModel(Vector<Core> MyCores) {
        this.MyCores = MyCores;
        if(this.MyCores==null){
            this.MyCores = new Vector<Core>();
        }
    }

    public Core getCore(int index){
        if(MyCores==null||MyCores.isEmpty()||index>=MyCores.size()){
            return null;
        }
        return MyCores.get(index);
    }

    @Override
    public int getRowCount() {
        return MyCores==null?0:MyCores.size();
    }

    @Override
    public int getColumnCount() {
        return columnNames.length;
    }

    @Override
    public String getColumnName(int column) {
        return columnNames[column];
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        if(column==0){
            return true;
        }else {
            return false;
        }
    }



    @Override
    public void setValueAt(Object aValue, int row, int column) {
        if(column==0 && aValue!=null&& aValue instanceof Boolean){
            MyCores.get(row).setSelected((Boolean)aValue);
//            fireTableDataChanged();
        }
    }

    @Override
    public Object getValueAt(int row, int column) {
        if(MyCores.size()==0){
            return null;
        }
        if(MyCores.size()<=row){
            return null;
        }
        switch(column){
            case 0:
                return MyCores.get(row).getSelected();
            case 1:
                return MyCores.get(row).getCore_mol();
            case 2:
                return MyCores.get(row).getName();
            case 3:
                return MyCores.get(row).getDescription();
            case 4:
                return MyCores.get(row).getChemist();
            case 5:
                return MyCores.get(row).getDate();
            case 6:
                return MyCores.get(row).getSource();
            case 7:
                return MyCores.get(row).getComment();
            case 8:
                return MyCores.get(row).getStatus();
            default:
                return null;
        }
    }

    @Override
    public Class<?> getColumnClass(int column) {
        switch(column){
            case 0:
                return Boolean.class;
            case 1:
                return PropertyMolecule.class;
            case 2:
                return String.class;
            case 3:
                return String.class;
            case 4:
                return String.class;
            case 5:
                return Date.class;
            case 6:
                return String.class;
            case 7:
                return String.class;
            case 8:
                return String.class;
            default:
                return Object.class;
        }
    }



}
