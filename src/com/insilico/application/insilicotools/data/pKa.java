package com.insilico.application.insilicotools.data;

/**
 * Created by jfeng1 on 10/5/15.
 */
public class pKa implements Comparable<pKa>{
    int atom_id;
    PropertyMolecule parentMol;
    float value;
    pKa_type type;

    public pKa(int atom_id, PropertyMolecule parentMol, float value, pKa_type type) {
        this.atom_id = atom_id;
        this.parentMol = parentMol;
        this.value = value;
        this.type = type;
    }



    @Override
    public int compareTo(pKa o) {
        if(type!=o.type){
            return 0;
        }
        if(type==pKa_type.ACIDIC){
            return new Float(value).compareTo(new Float(o.value));
        }else{
            return new Float(o.value).compareTo(new Float(value));
        }
    }

    @Override
    public String toString() {
        return(String.format("%s %d %5.2f",type,atom_id,value));
    }

    public int getAtom_id() {
        return atom_id;
    }

    public PropertyMolecule getParentMol() {
        return parentMol;
    }

    public float getValue() {
        return value;
    }

    public pKa_type getType() {
        return type;
    }

    public boolean isAcidic(){
        return type==pKa_type.ACIDIC;
    }
}
