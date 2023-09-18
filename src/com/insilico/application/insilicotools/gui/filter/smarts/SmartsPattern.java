package com.insilico.application.insilicotools.gui.filter.smarts;

public class SmartsPattern {
    String name;
    String smarts;

    public SmartsPattern(String name, String smarts) {
        this.name = name;
        this.smarts = smarts;
    }

    public String getName() {
        return name;
    }

    public String getSmarts() {
        return smarts;
    }

    @Override
    public String toString() {
        return name == null ? "" : name;
    }

    public static SmartsPattern[] getPredefinedSmartsPattern() {
        return new SmartsPattern[]{
                new SmartsPattern("acid chlorides", "[#6]C(=O)Cl"),
                new SmartsPattern("sulfony chlorides","[#6]S(=O)(=O)Cl"),
                new SmartsPattern("aldehydes", "[#6][CH1](=O)"),
                new SmartsPattern("alkyl boronates", "[CX4,CX3]B([OH])[OH]"),
                new SmartsPattern("aromatic bromides", "[#6;a]Br"),
                new SmartsPattern("aromatic chlorides", "[#6;a]Cl"),
                new SmartsPattern("aromatic fluorides", "[#6;a]F"),
                new SmartsPattern("aromatic iodides", "[#6;a]I"),
                new SmartsPattern("aryl boronates", "[#6;a]B([OH])[OH]"),
                new SmartsPattern("carboxylates", "[#6]C(=O)[OD1]"),
                new SmartsPattern("primary alcohols", "[CH2][OH]"),
                new SmartsPattern("primary amines", "[CX4,c][NH2]"),
                new SmartsPattern("primary bromides", "[CH2]Br"),
                new SmartsPattern("primary chlorides", "[CH2]Cl"),
                new SmartsPattern("primary iodides", "[CH2]I"),
                new SmartsPattern("secondary alcohols", "[#6][CX4;H1][OH]"),
                new SmartsPattern("secondary amines", "[CX4,c][NH1][CX4,c]"),
                new SmartsPattern("secondary bromides", "[CX4;H1]Br"),
                new SmartsPattern("secondary chlorides", "[CX4;H1]Cl"),
                new SmartsPattern("secondary iodides", "[CX4;H1]I"),
                new SmartsPattern("tertiary alcohols", "[CX4,c;H0][OH]"),
                new SmartsPattern("alcohols", "[CX4,c][OH]"),
                new SmartsPattern("primary and secondary amines", "[NX3;H2,H1;!$(NC=O);!$(NS=O)]"),
                new SmartsPattern("boronic acids and esters", "[#6]B(O)O"),
                new SmartsPattern("primary and secondary bromides", "[CX4;H2,H1]Br"),
                new SmartsPattern("primary and secondary chlorides", "[CX4;H2,H1]Cl"),
                new SmartsPattern("primary and secondary iodides", "[CX4;H2,H1]I")

        };
    }
}
