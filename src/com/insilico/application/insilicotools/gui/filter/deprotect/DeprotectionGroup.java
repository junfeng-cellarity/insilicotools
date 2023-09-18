package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.google.common.base.Function;

import java.io.Serializable;

public class DeprotectionGroup implements Serializable {
    private final Integer id;
    private final String name;
    private final String smirks;

    public static Function<DeprotectionGroup, String> depGroupToSmirksFunction = new Function<DeprotectionGroup, String>(){
        @Override
        public String apply(DeprotectionGroup deprotectionGroup) {
            return deprotectionGroup.getSmirks();
        }
    };

    public static Function<DeprotectionGroup, String> depGroupToNameFunction = new Function<DeprotectionGroup, String>(){
        @Override
        public String apply(DeprotectionGroup deprotectionGroup) {
            return deprotectionGroup.getName();
        }
    };
    private static DeprotectionGroup[] deprotectionGroups = new DeprotectionGroup[]{
            new DeprotectionGroup(1,"Boc","[#6:2]-[#7:1]-[#6](=O)-[#8]C([#6H3])([#6H3])[#6H3]>>[#6:2]-[#7:1]"),
            new DeprotectionGroup(2,"n-Boc","[#6:2]:[#7:1]-[#6](=O)-[#8]C([#6H3])([#6H3])[#6H3]>>[#6:2]-[#7:1]"),
            new DeprotectionGroup(3,"Cbz","[#6:2]-[#7:1]-[#6](=O)-[#8]-[#6H2]-c1[cH1][cH1][cH1][cH1][cH1]1>>[#6:2]-[#7:1]"),
            new DeprotectionGroup(4,"Fmoc","[#6:2]-[#7:1]-[#6](=O)-[#8]-[#6]-[#6]-1-c2ccccc2-c2ccccc-12>>[#6:2]-[#7:1]"),
            new DeprotectionGroup(5,"Benzyl ester","[#6:1]-[#6:2](=[O:4])-[#8:3]-[#6]-c1ccccc1>>[#6:1]-[#6:2](-[#8:3])=[O:4]"),
            new DeprotectionGroup(6,"T-butyl easter","[#6:1]-[#6:2](=[O:4])-[#8:3]C([#6])([#6])[#6]>>[#6:1]-[#6:2](-[#8:3])=[O:4]"),
            new DeprotectionGroup(7, "Bn","[#6:2]-[#7:1]-[#6H2]-c1[cH1][cH1][cH1][cH1][cH1]1>>[#6:2]-[#7:1]")
    };

    public DeprotectionGroup(Integer id, String name, String smirks) {
        this.id = id;
        this.name = name;
        this.smirks = smirks;
    }

    public static DeprotectionGroup[] getDeprotectionGroups() {
        return deprotectionGroups;
    }

    public String getName() {
        return name;
    }

    public String getSmirks() {
        return smirks;
    }

    public Integer getId() {
        return id;
    }
}
