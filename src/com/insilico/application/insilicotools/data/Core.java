package com.insilico.application.insilicotools.data;

import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;

import java.sql.SQLException;
import java.util.Date;

/**
 * Created by jfeng1 on 2/10/17.
 */
public class Core {
    int id;
    String core_mol;
    String name;
    String core_smiles;
    String description;
    byte[] cdx;
    String chemist;
    Date date;
    Boolean selected = false;
    String source;
    String comment;
    int status_id;
    String status;

    public Core(int id, String name, String core_mol, String core_smiles, String description, byte[] cdx, String chemist, Date date, String source, String comment, int status_id) {
        this.id = id;
        this.core_mol = core_mol;
        this.name = name;
        this.core_smiles = core_smiles;
        this.description = description;
        this.cdx = cdx;
        this.chemist = chemist;
        this.date = date;
        this.selected = false;
        this.source = source;
        this.comment = comment;
        this.status_id = status_id;
    }

    public void setCore_mol(String core_mol) {
        this.core_mol = core_mol;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setCore_smiles(String core_smiles) {
        this.core_smiles = core_smiles;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public void setCdx(byte[] cdx) {
        this.cdx = cdx;
    }

    public void setChemist(String chemist) {
        this.chemist = chemist;
    }

    public void setDate(Date date) {
        this.date = date;
    }

    public void setSource(String source) {
        this.source = source;
    }

    public void setComment(String comment) {
        this.comment = comment;
    }

    public void setStatus_id(int status_id) {
        this.status_id = status_id;
        this.status = retrieveStatus();
    }


    public String getComment() {
        return comment;
    }

    public int getStatus_id() {
        return status_id;
    }

    private String retrieveStatus(){
        try {
            return FrontierDAO.getInstance().getCoreStatus(status_id);
        } catch (SQLException e) {
            e.printStackTrace();
            return "";
        }
    }

    public String getStatus() {
        if(status==null){
            status = retrieveStatus();
        }
        return status;
    }

    public String getName() {
        return name;
    }

    public String getSource() {
        return source;
    }

    public Boolean getSelected() {
        return selected;
    }

    public void setSelected(Boolean selected) {
        this.selected = selected;
    }

    @Override
    public int hashCode() {
        return id;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Core that = (Core) o;
        return id == that.id;
    }

    public int getId() {
        return id;
    }

    public PropertyMolecule getCore_mol() {
        OEGraphMol mol = ChemFunc.getMolFromMolString(core_mol, OEFormat.SDF);
        return new PropertyMolecule(mol);
    }

    public String getCore_smiles() {
        return core_smiles;
    }

    public String getDescription() {
        return description;
    }

    public byte[] getCdx() {
        return cdx;
    }

    public String getChemist() {
        return chemist;
    }

    public Date getDate() {
        return date;
    }


}
