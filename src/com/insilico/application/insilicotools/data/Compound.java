package com.insilico.application.insilicotools.data;

import com.insilico.application.insilicotools.database.FrontierDAO;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.google.common.hash.HashCode;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;

import java.sql.SQLException;
import java.util.Date;

/**
 * Created by jfeng1 on 2/10/17.
 */
public class Compound {
    int id;
    String molecule;
    String name;
    String smiles;
    byte[] cdx;
    Date date;
    Boolean selected = false;
    int status_id;
    Status status;
    int project_id;
    String project_name;
    int chemist_id;
    int rank;
    boolean isChanged = false;
    boolean rankingChanged = false;
    public static String[] PROPERTIES_TO_USE = {"CLogP","molecular weight","2d PSA","CY-NUMBER"};
    public static Class[] PROPERTY_CLASS = {Double.class,Double.class,Double.class,String.class};
    PropertyMolecule propertyMolecule;
    boolean hasProperty = false;


    public boolean hasProperty() {
        return hasProperty;
    }

    public Compound(int id, String name, String mol, String smiles, byte[] cdx, Date date, int project_id, int status_id, int rank, int chemist_id) {
        this.id = id;
        this.molecule = mol;
        this.name = name;
        this.smiles = smiles;
        this.cdx = cdx;
        this.date = date;
        this.selected = false;
        this.status_id = status_id;
        this.project_id = project_id;
        this.rank = rank;
        this.chemist_id = chemist_id;
        updatePropertyMol();
        hasProperty = false;
    }

    public Compound(int id, String name, String mol, String smiles, byte[] cdx, Date date, int project_id, int status_id, int rank, int chemist_id, double clogp, double mw, double psa) {
        this.id = id;
        this.molecule = mol;
        this.name = name;
        this.smiles = smiles;
        this.cdx = cdx;
        this.date = date;
        this.selected = false;
        this.status_id = status_id;
        this.project_id = project_id;
        this.rank = rank;
        this.chemist_id = chemist_id;
        if(clogp==0.0&&mw==0.0&&psa==0.0){
            hasProperty = false;
            updatePropertyMol();
        }else {
            hasProperty = true;
            updatePropertyMol(clogp, mw, psa);
        }
    }

    void updatePropertyMol(double clogp, double mw, double psa) {
        OEGraphMol oemol = ChemFunc.getMolFromMolString(molecule, OEFormat.SDF);
        propertyMolecule = new PropertyMolecule(oemol);
        propertyMolecule.setName(this.name);
        propertyMolecule.addProperty("CLogP",new Double(clogp).toString());
        propertyMolecule.addProperty("molecular weight",new Double(mw).toString());
        propertyMolecule.addProperty("2d PSA", new Double(psa).toString());
    }

    public String getChemist(){
        if(chemist_id!=-1){
            return FrontierDAO.getInstance().getAssignedChemist(chemist_id).chemist_name;
        }else{
            return null;
        }
    }

    void updatePropertyMol() {
        OEGraphMol oemol = ChemFunc.getMolFromMolString(molecule, OEFormat.SDF);
        propertyMolecule = new PropertyMolecule(oemol);
        propertyMolecule.setName(this.name);
    }


    public Chemist getAssignedChemist(){
        if(chemist_id==0){
            return null;
        }
        return FrontierDAO.getInstance().getAssignedChemist(chemist_id);
    }

    public void setAssignedChemist(int chemist_id){
        this.chemist_id = chemist_id;
        isChanged = true;
    }

    public void setMol(String mol) {
        this.molecule = mol;
        updatePropertyMol();
        isChanged = true;
    }

    public boolean isChanged() {
        return isChanged;
    }

    public void setName(String name) {
        this.name = name;
        isChanged = true;
    }

    public boolean isAdded(){
        if(this.id==-1){
            return true;
        }
        return false;
    }


    public void setCdx(byte[] cdx) {
        this.cdx = cdx;
        isChanged = true;
    }

    public void setChemistId(int chemist_id) {
        this.chemist_id = chemist_id;
        isChanged = true;
    }

    public void setDate(Date date) {
        this.date = date;
        isChanged = true;
    }

    public void setStatus_id(int status_id) {
        this.status = null;
        this.status_id = status_id;
        this.status = getStatus();
        isChanged = true;
    }


    public int getStatus_id() {
        return status_id;
    }

    public Status getStatus() {
        if(status==null){
            status = FrontierDAO.getInstance().getStatus(status_id);
        }
        return status;
    }

    public String getName() {
        return name;
    }

    public Boolean getSelected() {
        return selected;
    }

    public void setSelected(Boolean selected) {
        this.selected = selected;
    }

    @Override
    public int hashCode() {
        if(id==-1){
            return HashCode.fromString(smiles).hashCode();
        }
        return id;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Compound that = (Compound) o;
        if(this.id==-1||that.id==-1){
            return false;
        }
        return id == that.id;
    }

    public int getId() {
        return id;
    }

    public PropertyMolecule getPropertyMol() {
        return propertyMolecule;
    }

    public String getSmiles() {
        return smiles;
    }

    public byte[] getCdx() {
        return cdx;
    }

    public int getChemistId() {
        return chemist_id;
    }

    public Date getDate() {
        return date;
    }

    public String getProject() {
        if(project_name==null){
            Project project = FrontierDAO.getInstance().getMyProject(project_id);
            if(project!=null) {
                project_name = project.getProject_name();
            }
        }
        return project_name;
    }

    public int getProject_id() {
        return project_id;
    }

    public int getRank() {
        return rank;
    }


    public void setRank(int rank) {
        this.rank = rank;
        isChanged = true;
        rankingChanged = true;
    }

    public void setChanged(boolean changed) {
        isChanged = changed;
    }

    public void setRankingChanged(boolean changed){
        rankingChanged = false;
    }

    public boolean isRankingChanged() {
        return rankingChanged;
    }

    public int save() throws SQLException {
        if(id==-1) {
                id = FrontierDAO.getInstance().insertMyCompound(this);
        }else{
            FrontierDAO.getInstance().updateMyCompound(this);
        }
        hasProperty = true;
        setChanged(false);
        setRankingChanged(false);
        return id;
    }
}
