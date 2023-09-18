package com.insilico.application.insilicotools.gui.widget;

import java.util.Vector;

public class MorphingFragment {
    String smiles;
    String name;
    String molFile;

    public MorphingFragment(String smiles, String name) throws Exception {
        this.smiles = smiles;
        this.name = name;
    }

    public MorphingFragment(String smiles, String name, String molFile) {
        this.smiles = smiles;
        this.name = name;
        this.molFile = molFile;
    }

    @Override
    public String toString() {
        return name;
    }

    public String getSmiles() {
        return smiles;
    }

    public String getName() {
        return name;
    }

    public String getMolFile() {
        return molFile;
    }

    public void setSmiles(String smiles) {
        this.smiles = smiles;
    }

    public void setName(String name) {
        this.name = name;
    }

    public void setMolFile(String molFile) {
        this.molFile = molFile;
    }

    /*
    public String getMolFile() {
        return molFile;
    }
    */

    public static Vector<MorphingFragment> getCommonGroups() {
        try {
            Vector<MorphingFragment> groups = new Vector<MorphingFragment>();
            groups.add(new MorphingFragment("C", "-CH3", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("Cl", "-Cl", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("C(F)(F)F", "-CF3", "\n  Marvin  11030916113D          \n" +
                    "\n" +
                    "  4  3  0  0  0  0            999 V2000\n" +
                    "   -0.0169    1.3948    0.0097 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
                    "    0.0021   -0.0041    0.0020 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    1.2956    1.8790   -0.0002 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "   -0.6935    1.8582   -1.1237 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "  1  2  1  0  0  0  0\n" +
                    "  1  3  1  0  0  0  0\n" +
                    "  1  4  1  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("F", "-F", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("O", "-OH", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("N", "-NH2", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("S(=O)(=O)N", "-SO2NH2", "\n" +
                    "  Marvin  11030916163D          \n" +
                    "\n" +
                    "  4  3  0  0  0  0            999 V2000\n" +
                    "   -0.0172    1.4168    0.0098 S   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    0.0021   -0.0041    0.0020 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    1.1617    2.2101    0.0028 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "   -0.9002    1.8797   -1.3124 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "  1  2  2  0  0  0  0\n" +
                    "  1  3  2  0  0  0  0\n" +
                    "  1  4  1  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("C(C)C", "-iPr", "\n" +
                    "  Marvin  11030916173D          \n" +
                    "\n" +
                    "  3  2  0  0  0  0            999 V2000\n" +
                    "   -0.0187    1.5258    0.0104 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    0.0021   -0.0041    0.0020 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    1.4167    2.0553   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "  1  2  1  0  0  0  0\n" +
                    "  1  3  1  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("C(C)(C)C", "-tBu", "\n" +
                    "  Marvin  11030916183D          \n" +
                    "\n" +
                    "  4  3  0  0  0  0            999 V2000\n" +
                    "   -0.0187    1.5258    0.0104 C   0  0  1  0  0  0  0  0  0  0  0  0\n" +
                    "    0.0021   -0.0041    0.0020 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    1.4167    2.0553   -0.0004 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "   -0.7587    2.0326   -1.2291 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "  1  2  1  0  0  0  0\n" +
                    "  1  3  1  0  0  0  0\n" +
                    "  1  4  1  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("OP(=O)(=O)O", "-OPO3", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  5  4  0  0  0  0            999 V2000\n" +
                    "   -0.0198    1.6058    0.0109 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    0.0021   -0.0041    0.0020 P   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    0.6948   -0.4814   -1.2156 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "    0.7179   -0.4943    1.2010 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "   -1.5084   -0.5613    0.0134 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "  1  2  1  0  0  0  0\n" +
                    "  2  3  2  0  0  0  0\n" +
                    "  2  4  2  0  0  0  0\n" +
                    "  2  5  1  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("I", "-I", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            groups.add(new MorphingFragment("Br", "-Br", "\n" +
                    "  Marvin  11040908563D          \n" +
                    "\n" +
                    "  1  0  0  0  0  0            999 V2000\n" +
                    "    0.0000    0.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n" +
                    "M  END"));
            return groups;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;

    }
}
