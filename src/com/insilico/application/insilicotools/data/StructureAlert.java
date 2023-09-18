package com.insilico.application.insilicotools.data;

import com.google.common.hash.HashCode;
import openeye.oechem.OEGraphMol;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;

/**
 * Created by jfeng1 on 10/9/15.
 */
public class StructureAlert {
    String smarts;
    String name;
    String description;
    int count = 0;
    OESubSearch subSearch;

    public StructureAlert(String smarts,String name, String description){
        this(smarts,name,description,0);
    }

    public StructureAlert(String smarts, String name, String description, int count) {
        this.smarts = smarts;
        this.name = name;
        this.description = description;
        subSearch = new OESubSearch(this.smarts);
        this.count = count;
        if(count>0){
            this.name = this.name+">"+count;
        }
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public OESubSearch getSubSearch() {
        return subSearch;
    }

    @Override
    public String toString() {
        return name + ":" + description;
    }

    @Override
    public boolean equals(Object obj) {
        if(obj instanceof StructureAlert){
            StructureAlert b = (StructureAlert)obj;
            return b.smarts.equals(smarts);
        }
        return false;
    }

    @Override
    public int hashCode() {
        return HashCode.fromString(smarts).hashCode();
    }

    public static void main(String[] args) {
        OESubSearch ss = new OESubSearch("c1ccccc1");
        OEGraphMol mol = new OEGraphMol();
        oechem.OESmilesToMol(mol,"c1ccccc1CC");
        oechem.OEPrepareSearch(mol,ss);
        ss.Match(mol);
    }
}
