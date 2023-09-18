package com.insilico.application.insilicotools.gui.filter.alert;

import com.insilico.application.insilicotools.InSilicoToolOptions;
import openeye.oechem.*;

import java.util.Vector;

public class AlertRule {
    private boolean selected = true;
    private String smarts;
    private int max;
    private String description;
    private int index;
    private OESubSearch subSearch;

    public AlertRule(String smarts, String description, int max, int index) {
        this.smarts = smarts;
        this.max = max;
        this.description = description;
        this.index = index;
        this.subSearch = new OESubSearch();
        this.subSearch.Init(smarts);
    }

    public boolean isPassed(OEGraphMol oemol){
        if(this.max>0) {
            return this.getMatchCount(oemol) <= this.max;
        }else{
            return !this.isSingleMatched(oemol);
        }
    }

    private boolean isSingleMatched(OEGraphMol oemol){
        if(oemol!=null){
            boolean hydrogenAdded = false;
            if(!oechem.OEHasExplicitHydrogens(oemol)) {
                hydrogenAdded = oechem.OEAddExplicitHydrogens(oemol);
            }
            oechem.OEPrepareSearch(oemol,this.subSearch);
            boolean singleMatch = this.subSearch.SingleMatch(oemol);
            if(hydrogenAdded) {
                oechem.OESuppressHydrogens(oemol);
            }
            return singleMatch;
        }
        return false;
    }

    private int getMatchCount(OEGraphMol oemol){
        if(oemol!=null){
            boolean hydrogenAdded = false;
            if(!oechem.OEHasExplicitHydrogens(oemol)) {
                hydrogenAdded = oechem.OEAddExplicitHydrogens(oemol);
            }
            oechem.OEPrepareSearch(oemol,this.subSearch);
            int matchCount = 0;
            for(OEMatchBase ma:this.subSearch.Match(oemol)){
                if(ma.IsValid()){
                    matchCount += 1;
                }
            }
            if(hydrogenAdded){
                oechem.OESuppressHydrogens(oemol);
            }
            return matchCount;
        }
        return 0;
    }

    public void setMax(int max) {
        this.max = max;
    }

    public boolean isSelected() {
        return selected;
    }

    public void setSelected(final boolean selected) {
        this.selected = selected;
    }

    public String getSmarts() {
        return smarts;
    }

    public int getMax() {
        return max;
    }

    public String getDescription() {
        return description;
    }

    public int getIndex() {
        return index;
    }

    @Override
    public String toString() {
        return description;
    }

    public static Vector<AlertRule> generatePainsRules(){
        Vector<AlertRule> alertRules = new Vector<AlertRule>();
        String pains_smarts = InSilicoToolOptions.pains_smarts;
        String[] smarts = pains_smarts.split("\n");
        int idx = 0;
        for(String s :smarts){
            if(s.startsWith("#")){
                continue;
            }
            String[] p = s.split("\t");
            if(p.length==2){
                alertRules.add(new AlertRule(p[1],p[0], 0, idx ++));
            }
        }

        return alertRules;
    }


    public static AlertRule[] LIBRARY_ALERT_RULES = new AlertRule[]{
            new AlertRule("[F]","Fluoro",6,0),
            new AlertRule("[Br]","Bromo",1,1),
            new AlertRule("[N+](=[OX1])[OX1]","Nitro",0,2),
            new AlertRule("[CX3;!r]=[CX3;!r]","C=C",0,3),
            new AlertRule("[CX3](=[OX1])[CX3](=[OX1])","Oxaldehyde",0,4),
            new AlertRule("[CX3](=[OX1])[CX4]([FX1])([FX1])[FX1]","TFA",0,5),
            new AlertRule("[Cl,Br,I]","Cl,Br,I",3,6),
            new AlertRule("[I]","I",0,7),
            new AlertRule("O-O","O-O",0,8),
            new AlertRule("S-S","S-S",0,9),
            new AlertRule("[N]-[N]","N-N",0,10),
            new AlertRule("[N]=[N]","N=N",0,11),
            new AlertRule("[N][N;H2]","NNH2",0,12),
            new AlertRule("[CX4;!$(C~[N,O,S]);!r][CX4;!$(C~[N,O,S]);!r][CX4;!$(C~[N,O,S]);!r][CX4;!$(C~[N,O,S]);!r]","Long Carbon Chain",0,13),
            new AlertRule("[!#1;!#6;!#7;!#8;!#9;!#35;!#17;!#53;!#16]","Uncommon Elements",0,14),
            new AlertRule("[N][S;!$(S(=O)(=O)N)]","N-S",0,15),
            new AlertRule("S(~O)(~O)(~O)","SO3",0,16),
            new AlertRule("[N,S][Cl,F,I]","NX or SX",0,17),
            new AlertRule("a-[O;H1]","phenol",1,18),
            new AlertRule("a-[NX3;H2]","Aniline",1,19),
            new AlertRule("[A;!r]=[A;!r]-[A;!r]=[A;!r]","Michael Receptor",0,20),
            new AlertRule("a-[SX2;H1]","ArSH",0,21),
            new AlertRule("[CX4]-[Cl,Br,I]","Alkyl Halide",0,22),
            new AlertRule("C#C-C#C","diacetylene",0,23),
            new AlertRule("[N;!r]=[C;!r]-Cl","N=C-Cl",0,24),
            new AlertRule("*=C(C-O)C=O","hydroxypyruvaldehyde",0,25),
            new AlertRule("[N]#[N+]","diazynium",0,26),
            new AlertRule("[N,S,O]=C=N","N=C=N,O=C=N,S=C=N",0,27),
            new AlertRule("S-C#N","thiocyanate",0,28),
            new AlertRule("C=N=N","C=N=N",0,29),
            new AlertRule("C-N=O","C-N=O",0,30),
            new AlertRule("N~N~N","N3",0,31),
            new AlertRule("[#6]-[CX3;H1]=O","Aldehyde",0,32),
            new AlertRule("[C;!r]-[C;!r]=[N;!r]-[C;!r]","ethylidene(methyl)amine",0,33),
            new AlertRule("C(C#N)C#N","malononitrile",0,34),
            new AlertRule("[S;!r]=C","S=C",0,35),
            new AlertRule("C(=O)[F,Cl,Br,I]","Acetyl Halide",0,36),
            new AlertRule("C(=O)C#N","cyanoketone",0,37),
            new AlertRule("O=C-O-C=O","formic anhydride",0,38),
            new AlertRule("a-O-C(=O)A","ArylEster",0,39),
            new AlertRule("N=C-Cl","N=C-Cl",0,40),
            new AlertRule("[OD1&!R]=C1C(=[OD1&!R])C=CC=C1","1,2-quinone",0,41),
            new AlertRule("O=[#6]1[#6]:,=[#6][#6](=O)[#6]:,=[#6]1","1,4-quinone",0,42),
            new AlertRule("[CX4;r][OX2;H1].[CX4;r][OX2;H1].[CX4;r][OX2;H1].[CX4;r][OX2;H1]","Sugar Like",0,43),
            new AlertRule("[#7,#8,#16;r3]","epoxide like",0,44),
            new AlertRule("[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N[a])]","Amine",1,45),
            new AlertRule("*(=O)-[O;!r]-*","Ester or Acid", 0, 46),
            new AlertRule("c1cc[#8,#16;a]c1","Furan Or Thiophene", 0, 47),
    };

    public static void main(String[] args) {
        if(args.length==2){
            String input_file = args[0];
            String output_file = args[1];
            oemolistream ifs = new oemolistream();
            ifs.open(input_file);
            oemolostream ofs = new oemolostream();
            ofs.open(output_file);
            OEGraphMol mol = new OEGraphMol();
            while(oechem.OEReadMolecule(ifs,mol)){
                boolean passed = true;
                for(AlertRule rule : LIBRARY_ALERT_RULES){
                    if(!rule.isPassed(mol)){
                        passed = false;
                        break;
                    }
                }
                if(passed){
                    oechem.OEWriteMolecule(ofs, mol);
                }
            }
            ifs.close();
            ofs.close();

        }else{
            System.out.println("missing input file.");
        }
    }
}