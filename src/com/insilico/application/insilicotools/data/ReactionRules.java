package com.insilico.application.insilicotools.data;

import chemaxon.struc.PeriodicSystem;

/**
 * Created by jfeng1 on 3/4/16.
 */
public class ReactionRules {
    public static ReactionRules ANY = new ReactionRules(-1,"ANY","");
    public static ReactionRules PRIMARY_AMINE_NITROGEN = new ReactionRules(PeriodicSystem.N, "Primary Amine","[#7;X3;H2;!$(NC=O)]");
    public static ReactionRules SECONDARY_AMINE_NITROGEN = new ReactionRules(PeriodicSystem.N, "Secondary Amine","[#7;X3;H1;!$(NC=O)] ");
    public static ReactionRules PRIMARY_OR_SECONDARY_AMINE_NITROGEN = new ReactionRules(PeriodicSystem.N, "Primary or Secondary Amine","[#7;X3;H2,H1;!$(NC=O)]");

    public static ReactionRules OXYGEN_SP3_ONE_HYDROGEN = new ReactionRules(PeriodicSystem.O, "Sp3, one hydrogen","[#8;X2;H1]");
    public static ReactionRules OXYGEN_SP3_NO_HYDROGEN = new ReactionRules(PeriodicSystem.O, "Sp3, no hydrogen","[#8;X2;H0]");
    public static ReactionRules OXYGEN_SP2 = new ReactionRules(PeriodicSystem.O, "Sp2","[#8;X1]");

    public static ReactionRules CARBON_AROMATIC = new ReactionRules(PeriodicSystem.C,"Aromatic Carbon","[#6;a]");
    public static ReactionRules CARBON_ALIPHATIC = new ReactionRules(PeriodicSystem.C,"Aliphatic Carbon","[#6;A]");


    public static ReactionRules[] getNitrogenRules(){
        return new ReactionRules[]{ANY, PRIMARY_OR_SECONDARY_AMINE_NITROGEN,PRIMARY_AMINE_NITROGEN,SECONDARY_AMINE_NITROGEN};
    }

    public static ReactionRules[] getOxygenRules(){
        return new ReactionRules[]{ANY, OXYGEN_SP2,OXYGEN_SP3_ONE_HYDROGEN,OXYGEN_SP3_NO_HYDROGEN};
    }

    public static ReactionRules[] getCarbonRules(){
        return new ReactionRules[]{ANY, CARBON_ALIPHATIC, CARBON_AROMATIC};
    }

    int elementno;
    String name;
    String smarts;

    public ReactionRules(int elementno, String name, String smarts) {
        this.elementno = elementno;
        this.name = name;
        this.smarts = smarts;
    }

    @Override
    public String toString() {
        return name;
    }

    public boolean isAny(){
        return name.equals("ANY");
    }

    public int getElementno() {
        return elementno;
    }

    public String getName() {
        return name;
    }

    public String getSmarts() {
        return smarts;
    }
}
