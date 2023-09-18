package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.*;

import javax.swing.table.DefaultTableModel;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;

/**
 * Created by jfeng1 on 2/13/17.
 */
public class MyCompoundTableModel extends DefaultTableModel{
    Vector<Compound> MyCompounds;
    String[] columnNames = new String[]{"","Structure","Name", "Chemist","Date", "Status", "Project", "Rank", "Assign To"};

    public MyCompoundTableModel(Vector<Compound> MyCompounds) {
        this.MyCompounds = MyCompounds;
        if(this.MyCompounds==null){
            this.MyCompounds = new Vector<Compound>();
        }
    }

    public Compound getCompoundByIndex(int index){
        if(MyCompounds==null||MyCompounds.isEmpty()||index>=MyCompounds.size()){
            return null;
        }
        return MyCompounds.get(index);
    }

    @Override
    public int getRowCount() {
        return MyCompounds==null?0:MyCompounds.size();
    }

    @Override
    public int getColumnCount() {
        return columnNames.length+ Compound.PROPERTIES_TO_USE.length;
    }

    @Override
    public String getColumnName(int column) {
        if(column<columnNames.length) {
            return columnNames[column];
        }else{
            return Compound.PROPERTIES_TO_USE[column-columnNames.length];
        }
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        if(column==0){
            return true;
        }else if(column==5){ //status
            return true;
        }else if(column==7){
            return true;
        }else if(column==8){
            return true;
        }
        else {
            return false;
        }
    }

    @Override
    public void setValueAt(Object aValue, int row, int column) {
        if(column==0 && aValue!=null&& aValue instanceof Boolean){
            MyCompounds.get(row).setSelected((Boolean)aValue);
//            fireTableDataChanged();
        }else if(column==5){
            MyCompounds.get(row).setStatus_id(((Status)aValue).getStatus_id());
        }else if(column==7){
            Integer rank = (Integer) aValue;
            MyCompounds.get(row).setRank(rank);
        }else if(column==8){
            Chemist chemist = (Chemist) aValue;
            if(chemist!=null) {
                MyCompounds.get(row).setAssignedChemist(chemist.getChemist_id());
            }
        }
    }

    @Override
    public Object getValueAt(int row, int column) {
        if(MyCompounds.size()==0){
            return null;
        }
        if(MyCompounds.size()<=row){
            return null;
        }
        switch(column){
            case 0:
                return MyCompounds.get(row).getSelected();
            case 1:
                return MyCompounds.get(row).getPropertyMol();
            case 2:
                return MyCompounds.get(row).getName();
            case 3:
                return MyCompounds.get(row).getChemist();
            case 4:
                return MyCompounds.get(row).getDate();
            case 5:
                return MyCompounds.get(row).getStatus();
            case 6:
                return MyCompounds.get(row).getProject();
            case 7:
                return MyCompounds.get(row).getRank();
            case 8:
                return MyCompounds.get(row).getAssignedChemist();
            default:
                PropertyMolecule pmol = MyCompounds.get(row).getPropertyMol();
                MolProperty property = pmol.getProperty(Compound.PROPERTIES_TO_USE[column - columnNames.length]);
                if(property!=null) {
                    if(Compound.PROPERTY_CLASS[column-columnNames.length].equals(Double.class)) {
                        return property.getValue();
                    }else{
                        return property.getProperty();
                    }
                }else{
                    return null;
                }
        }
    }

    public String getValueAsString(int row, int column) {
        if(MyCompounds.size()==0){
            return null;
        }
        if(MyCompounds.size()<=row){
            return null;
        }
        switch(column){
            case 0:
                return MyCompounds.get(row).getSelected().toString();
            case 1:
                return MyCompounds.get(row).getPropertyMol().getSmiles();
            case 2:
                return MyCompounds.get(row).getName();
            case 3:
                return MyCompounds.get(row).getChemist();
            case 4:
                SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy");
                return sdf.format(MyCompounds.get(row).getDate());
            case 5:
                return MyCompounds.get(row).getStatus().getStatus_name();
            case 6:
                return MyCompounds.get(row).getProject();
            case 7:
                return new Integer(MyCompounds.get(row).getRank()).toString();
            case 8:
                return MyCompounds.get(row).getAssignedChemist().getChemist_name();
            default:
                PropertyMolecule pmol = MyCompounds.get(row).getPropertyMol();
                MolProperty property = pmol.getProperty(Compound.PROPERTIES_TO_USE[column - columnNames.length]);
                if(property!=null) {
                    return property.toString();
                }else{
                    return null;
                }
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
                return Date.class;
            case 5:
                return Status.class;
            case 6:
                return String.class;
            case 7:
                return Integer.class;
            case 8:
                return Chemist.class;
            default:
                return Compound.PROPERTY_CLASS[column-columnNames.length];
        }
    }



}
