package com.insilico.application.insilicotools.gui.lims;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import javax.swing.table.DefaultTableModel;
import java.util.Vector;

public class LimsMolTableModel extends DefaultTableModel {
    Vector<LimsMolecule> molecules;
    public LimsMolTableModel(Vector<LimsMolecule> mols) {
        if(mols!=null){
            molecules=mols;
        }else{
            molecules = new Vector<>();
        }
    }

    @Override
    public boolean isCellEditable(int row, int column) {
        if(column==1){
            return false;
        }
        if(column==2){
            return false;
        }
        return true;
    }



    @Override
    public void setValueAt(Object aValue, int row, int column) {
        OEGraphMol mol = molecules.get(row).getPropertyMolecule().getMol();
        switch (column){
            case 0:
                molecules.get(row).setSelected((Boolean)aValue);
                break;
            case LimsMolecule.COMPOUND_ID_IDX:
                molecules.get(row).setCompound_id((String)aValue);
                break;
            case LimsMolecule.NAME_IDX:
                molecules.get(row).setName((String)aValue);
                break;
            case LimsMolecule.THEORY_MASS_IDX:
                molecules.get(row).setTheory_mass((String)aValue);
                break;
            case LimsMolecule.PROJECT_IDX:
                molecules.get(row).setProject((String) aValue);
                break;
            case LimsMolecule.SCIENTIST_IDX:
                molecules.get(row).setScientist((String)aValue);
                break;
            default:
                if(oechem.OEHasSDData(mol,LimsMolecule.required_tags[column])) {
                    oechem.OESetSDData(mol, LimsMolecule.required_tags[column], (String) aValue);
                }else{
                    oechem.OEAddSDData(mol, LimsMolecule.required_tags[column], (String) aValue);
                }
                break;
        }
        molecules.get(row).save();
    }

    @Override
    public String getColumnName(int column) {
        return LimsMolecule.required_tags[column];
    }

    @Override
    public int getColumnCount() {
        return LimsMolecule.required_tags.length;
    }

    @Override
    public int getRowCount() {
        return molecules==null?0:molecules.size();
    }

    @Override
    public Class<?> getColumnClass(int column) {
        switch (column){
            case 0:
                return Boolean.class;
            case 1:
                return PropertyMolecule.class;
            case LimsMolecule.PROJECT_IDX:
                return String.class;
            default:
                return String.class;
        }
    }

    @Override
    public Object getValueAt(int row, int column) {
        if(row>=molecules.size()){
            return "";
        }
        OEGraphMol mol = molecules.get(row).getPropertyMolecule().getMol();
        switch (column){
            case 0:
                return molecules.get(row).getSelected();
            case 1:
                return molecules.get(row).getPropertyMolecule();
            case LimsMolecule.COMPOUND_ID_IDX:
                return molecules.get(row).getCompound_id();
            case LimsMolecule.NAME_IDX:
                return molecules.get(row).getName();
            case LimsMolecule.THEORY_MASS_IDX:
                return molecules.get(row).getTheory_mass();
            case LimsMolecule.PROJECT_IDX:
                return molecules.get(row).getProject();
            case LimsMolecule.SCIENTIST_IDX:
                return molecules.get(row).getScientist();
            default:
                if(oechem.OEHasSDData(mol,LimsMolecule.required_tags[column])){
                    return oechem.OEGetSDData(mol,LimsMolecule.required_tags[column]);
                }else{
                    return "";
                }
        }
    }
}
