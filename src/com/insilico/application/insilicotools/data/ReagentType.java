package com.insilico.application.insilicotools.data;

/**
 * Created by jfeng1 on 8/3/16.
 */
public class ReagentType {
    public static final String PRIMARY_ALKYL_AMINE = "[#6]-[#7;X3;H2;!$(NC=O);!$(Nc)]";
    public static final String SECONDARY_ALKYL_AMINE = "[#6]-[#7;X3;H1;!$(NC=O);!$(Nc)]";
    public static final String PRIMARY_SECONDARY_ALKYL_AMINE = "[#6]-[#7;X3;H2,H1;!$(NC=O);!$(Nc)]";
    public static final String PRIMARY_AMINE = "[#6]-[#7;X3;H2;!$(NC=O)]";
    public static final String SECONDARY_AMINE = "[#6]-[#7;X3;H1;!$(NC=O)]";
    public static final String PRIMARY_SECONDARY_AMINE = "[#6]-[#7;X3;H2,H1;!$(NC=O)]";
    public static final String ALDEHYDE = "[#6]-[#6H1]=O";
    public static final String SULFONYL_CHLORIDE = "[#6][S](Cl)(=[O])=[O]";
    public static final String CHLOROFORMATE = "[#6]-[#8]-[#6](Cl)=[O]";
    public static final String CARBOXYLIC_ACID = "[#6]-[#6](-[#8H1])=[O]";
    public static final String ARYL_HALIDE = "[#6;a]-[#17,#35,#53]";
    public static final String ALKYL_HALIDE = "[#6;A]-[#17,#35,#53]";
    public static final String HALIDE = "[#6]-[#17,#35,#53]";
    public static final String BORONIC_ACIDS_AND_ESTERS = "[#6]B(O)O";

}
