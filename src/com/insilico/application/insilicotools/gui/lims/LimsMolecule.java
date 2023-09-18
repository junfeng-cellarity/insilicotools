package com.insilico.application.insilicotools.gui.lims;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.database.LimsDAO;
import com.google.common.hash.HashCode;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.HashMap;

public class LimsMolecule implements Comparable<LimsMolecule>{
    int cmpd_id;
    String name;
    String compound_id;
    String scientist;
    Boolean selected = false;
    String project;
    String theory_mass;
    Timestamp timestamp;


    PropertyMolecule propertyMolecule;
    static final int COMPOUND_ID_IDX = 4;
    static final int NAME_IDX = 6;
    static final int THEORY_MASS_IDX = 10;
    static final int PROJECT_IDX = 16;
    static final int SCIENTIST_IDX = 17;
    static final int EM_IDX = 7;
    static final int USER_COMPOUND_ID_IDX = 21;
    static final int COMPATIBILITY_WITH_TFA_IDX = 18;
    static final int AMOUNT_SUBMITTED_IDX = 19;
    public static String[] required_tags = {"","Structure","Select","Pos.","Compound ID","Salt","Name","EM",
                                            "Mass","Succeed","Theory Mass","Vol.","Density","Purity","Moles",
                                            "Yield","Project","Scientist","Compatibility with TFA","amount submitted",
                                            "Note",	"User_CompoundID","DissolutionVol","InjectionNo","InjectionVol"
                                            };

    public LimsMolecule(PropertyMolecule propertyMolecule, int cmpd_id, String compound_id, String theory_mass, String project, String scientist, String compatibility_with_tfa, String amount_submitted, Long timestamp){
        this.propertyMolecule = propertyMolecule;
        if(timestamp==null) {
            this.timestamp = new Timestamp(0L);
        }else{
            this.timestamp = new Timestamp(timestamp);
        }
        if(propertyMolecule==null){
            this.propertyMolecule = new PropertyMolecule(new OEGraphMol());
        }
        HashMap<String,String> dataDict = new HashMap<>();
        this.compound_id = compound_id;
        if(compound_id==null){
            this.compound_id = "";
        }
        dataDict.put(required_tags[COMPOUND_ID_IDX],this.compound_id);
        dataDict.put(required_tags[USER_COMPOUND_ID_IDX],this.compound_id);
        if(compatibility_with_tfa!=null){
            dataDict.put(required_tags[COMPATIBILITY_WITH_TFA_IDX],compatibility_with_tfa);
        }
        if(amount_submitted!=null){
            dataDict.put(required_tags[AMOUNT_SUBMITTED_IDX],amount_submitted);
        }

        this.cmpd_id = cmpd_id;
        this.theory_mass = theory_mass;
        if(theory_mass==null){
            this.theory_mass = "";
        }
        dataDict.put(required_tags[THEORY_MASS_IDX],this.theory_mass);

        this.project = project;
        if(project == null){
            this.project = "";
        }
        dataDict.put(required_tags[PROJECT_IDX],this.project);

        this.scientist = scientist;
        if(scientist==null){
            this.scientist = "";
        }
        dataDict.put(required_tags[SCIENTIST_IDX],this.scientist);

        this.name = propertyMolecule.getName();
        OEGraphMol oemol = propertyMolecule.getMol();
        dataDict.put(required_tags[EM_IDX], ""+ oechem.OECalculateMolecularWeight(oemol,true));
        for(String tag:LimsMolecule.required_tags){
            if(tag.equals("")||tag.equals("Structure")){
                continue;
            }
            if(oechem.OEHasSDData(oemol,tag)&&!oechem.OEGetSDData(oemol,tag).isEmpty()){
                dataDict.put(tag,oechem.OEGetSDData(oemol,tag));
            }
            oechem.OEDeleteSDData(oemol,tag);
        }
        for(String tag:LimsMolecule.required_tags){
            if(tag.equals("")||tag.equals("Structure")){
                continue;
            }
            String data = "";
            if(dataDict.containsKey(tag)){
                data = dataDict.get(tag);
            }
            oechem.OEAddSDData(oemol,tag,data);
        }
    }

    public void generateCompoundId(int batch_id, int id){
        OEGraphMol oemol = propertyMolecule.getMol();
        if(batch_id==-1){
            this.compound_id = oechem.OEGetSDData(oemol,required_tags[USER_COMPOUND_ID_IDX]);
        }else {
            this.compound_id = String.format("900000-%d-P%d", batch_id, id);
        }
        oechem.OESetSDData(oemol,required_tags[COMPOUND_ID_IDX],this.compound_id);
    }


    public String getCompound_id() {
        return compound_id;
    }

    public String getTheory_mass() {
        return theory_mass;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setScientist(String scientist) {
        if(scientist==null){
            scientist = "";
        }
        this.scientist = scientist;
        oechem.OESetSDData(propertyMolecule.getMol(),required_tags[SCIENTIST_IDX],scientist);
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
        if(cmpd_id==-1){
            if(propertyMolecule!=null) {
                if(propertyMolecule.getSdfStr()!=null) {
                    return HashCode.fromString(propertyMolecule.getSdfStr()).hashCode();
                }else{
                    return HashCode.fromString("").hashCode();
                }
            }else{
                return HashCode.fromString("").hashCode();
            }
        }
        return cmpd_id;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        LimsMolecule that = (LimsMolecule) o;
        if(this.cmpd_id==-1||that.cmpd_id==-1){
            return false;
        }
        return cmpd_id == that.cmpd_id;
    }

    public int getCmpdId() {
        return cmpd_id;
    }

    public String getScientist() {
        return scientist;
    }

    public String getProject() {
        return project;
    }

    public int getCmpd_id() {
        return cmpd_id;
    }


    public PropertyMolecule getPropertyMolecule() {
        return propertyMolecule;
    }

    public void setCmpd_id(int cmpd_id) {
        this.cmpd_id = cmpd_id;
    }

    public void setTheory_mass(String theory_mass) {
        if(theory_mass==null){
            theory_mass = "";
        }
        this.theory_mass = theory_mass;
        oechem.OESetSDData(propertyMolecule.getMol(),required_tags[THEORY_MASS_IDX],this.theory_mass);
    }

    public void setCompound_id(String compound_id) {
        if(compound_id==null){
            compound_id = "";
        }
        this.compound_id = compound_id;
        oechem.OESetSDData(propertyMolecule.getMol(),required_tags[COMPOUND_ID_IDX],this.compound_id);
    }

    public void setProject(String project) {
        if(project == null){
            project = "";
        }
        this.project = project;
        oechem.OESetSDData(propertyMolecule.getMol(),required_tags[PROJECT_IDX],this.project);
    }

    public Long getTimestamp() {
        return timestamp.getTime();
    }

    public void save(){
        try {
            if(cmpd_id==-1) {
                cmpd_id = LimsDAO.getInstance().insertLimsCompound(this);
            }else{
                LimsDAO.getInstance().updateLimsCompound(this);
            }
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    int extractInt(String s) {
        String num = s.replaceAll("\\D", "");
        // return 0 if no digits found
        return num.isEmpty() ? 0 : Integer.parseInt(num);
    }

    private int compare_compound_id(String compound_id_1, String compound_id_2){
        if(compound_id_1==null||compound_id_2==null){
            return 0;
        }
        if(compound_id_1.isEmpty()||compound_id_2.isEmpty()){
            return 0;
        }
        String[] arg_1 = compound_id_1.split("-");
        String[] arg_2 = compound_id_2.split("-");
        if(arg_1.length!=3||arg_2.length!=3){
            return 0;
        }

        int c1_0 = extractInt(arg_1[0]);
        int c2_0 = extractInt(arg_2[0]);

        int c1_1 = extractInt(arg_1[1]);
        int c2_1 = extractInt(arg_2[1]);

        int c1_2 = extractInt(arg_1[2]);
        int c2_2 = extractInt(arg_2[2]);

        if(c1_0!=c2_0){
            return c1_0-c2_0;
        }
        if(c1_1!=c2_1){
            return c1_1-c2_1;
        }
        return c1_2-c2_2;
    }


    @Override
    public int compareTo(LimsMolecule other) {
        if(this.timestamp.equals(other.timestamp)){
            return compare_compound_id(this.compound_id,other.compound_id);
        }else{
            return this.timestamp.compareTo(other.timestamp);
        }
    }

//    public static int extractInt2(String s) {
//        String num = s.replaceAll("\\D", "");
//        // return 0 if no digits found
//        return num.isEmpty() ? 0 : Integer.parseInt(num);
//    }
//    public static int test_compare_compound_id(String compound_id_1, String compound_id_2){
//        if(compound_id_1==null||compound_id_2==null){
//            return 0;
//        }
//        if(compound_id_1.isEmpty()||compound_id_2.isEmpty()){
//            return 0;
//        }
//        String[] arg_1 = compound_id_1.split("-");
//        String[] arg_2 = compound_id_2.split("-");
//        if(arg_1.length!=3||arg_2.length!=3){
//            return 0;
//        }
//
//        int c1_0 = extractInt2(arg_1[0]);
//        int c2_0 = extractInt2(arg_2[0]);
//
//        int c1_1 = extractInt2(arg_1[1]);
//        int c2_1 = extractInt2(arg_2[1]);
//
//        int c1_2 = extractInt2(arg_1[2]);
//        int c2_2 = extractInt2(arg_2[2]);
//
//        if(c1_0!=c2_0){
//            return c1_0-c2_0;
//        }
//        if(c1_1!=c2_1){
//            return c1_1-c2_1;
//        }
//        return c1_2-c2_2;
//    }

//    public static void main(String[] args) {
//        String x1 = "100228-80-1";
//        String x2 = "100228-80-2";
//        System.out.println(LimsMolecule.test_compare_compound_id(x1,x2));
//    }
}
