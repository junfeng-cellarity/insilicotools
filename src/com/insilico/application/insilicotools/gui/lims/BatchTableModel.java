package com.insilico.application.insilicotools.gui.lims;

import com.insilico.application.insilicotools.database.LimsDAO;

import javax.swing.table.DefaultTableModel;
import java.util.Date;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created by jfeng1 on 3/23/17.
 */
public class BatchTableModel extends DefaultTableModel {
    Vector<Batch> batches;
    String[] colNames = {"Batch Name", "Chemist", "Date", "No. Mols"};
    HashMap<Integer,Integer> batchNumMolsDict = new HashMap<>();

    public BatchTableModel(Vector<Batch> batches) {
        this.batches = batches;
        if(batches==null){
            this.batches = new Vector<Batch>();
        }
    }

    public void clearCache(){
        batchNumMolsDict.clear();
    }

    @Override
    public int getRowCount() {
        return batches==null?0:batches.size();
    }

    @Override
    public int getColumnCount() {
        return colNames.length;
    }

    @Override
    public String getColumnName(int column) {
        return colNames[column];
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        return false;
    }

    @Override
    public Object getValueAt(int row, int column) {
        switch(column){
            case 0:
                return batches.get(row).getBatchName();
            case 1:
                return batches.get(row).getChemist();
            case 2:
                return batches.get(row).getDate();
            case 3:
                if(!batchNumMolsDict.containsKey(batches.get(row).batch_id)){
                    batchNumMolsDict.put(batches.get(row).batch_id,LimsDAO.getInstance().getNumMolsByBatch(batches.get(row).getBatch_id()));
                }
                return batchNumMolsDict.get(batches.get(row).batch_id);
            default:
                return "";
        }
    }

    @Override
    public Class<?> getColumnClass(int columnIndex) {
        switch(columnIndex){
            case 0:
                return String.class;
            case 1:
                return String.class;
            case 2:
                return Date.class;
            case 3:
                return Integer.class;
            default:
                return Object.class;
        }
    }


}
