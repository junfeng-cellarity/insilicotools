package com.insilico.application.insilicotools.data;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.*;
import org.slf4j.LoggerFactory;

import java.util.*;

/**
 * Created by jfeng1 on 9/8/15.
 */
public class PropertyMolecule implements Comparable{
    OEGraphMol mol;
    OEGraphMol mol3d = null;
    String name;
    LinkedHashMap<String,MolProperty> propertyMap;
    SVGString svgString;
    boolean isSelected;
    OEAtomBondSet structure_alert_set;
    int heavyAtomCount = 0;
    String smiles;
    String uniq_mol_id;
    boolean isMarked = false;
    String formula = "";

//    private final static String[] OE_PROPERTIES = {
//            "hydrogen-bond acceptors",
//            "hydrogen-bond donors",
//            "molecular weight",
//            "sum of formal charges",
//            "number of rings",
//            "rotatable bonds",
//            "XLogP",
//            "Solubility",
//            "2d PSA",
//            "Lipinski violations"
//    };

    public final static String OE_HBA = ChemFunc.OE_PROPERTIES[0];
    public final static String OE_HBD = ChemFunc.OE_PROPERTIES[1];
    public final static String OE_MW = ChemFunc.OE_PROPERTIES[2];
    public final static String OE_SUM_CHARGES = ChemFunc.OE_PROPERTIES[3];
    public final static String OE_NO_RINGS = ChemFunc.OE_PROPERTIES[4];
    public final static String OE_RTB = ChemFunc.OE_PROPERTIES[5];
    public final static String OE_XLOGP = ChemFunc.OE_PROPERTIES[6];
    public final static String OE_TPSA = ChemFunc.OE_PROPERTIES[7];
    boolean needUpdate = false;
    public final static int CHEMAXON_PKA = 1;
    public final static int MOKA_PKA = 2;
    private LinkedHashMap<Integer,pKa> chemaxon_pKa_dict;
    private LinkedHashMap<Integer,pKa> moka_pKa_dict;
    private volatile static int mol_id = 1;

    public LinkedHashMap<Integer, pKa> getpKaMap(int type) {
        if(type==CHEMAXON_PKA){
            return chemaxon_pKa_dict;
        }else{
            return moka_pKa_dict;
        }
    }

    public void addPKa(int atomId, pKa p, int type){
        if(type == CHEMAXON_PKA){
            chemaxon_pKa_dict.put(atomId,p);
        }else{
            moka_pKa_dict.put(atomId,p);
        }
        needUpdate = true;
    }

    public void addStructureAlerts(StructureAlert[] alerts){
        if(alerts!=null){
            StringBuilder sb = new StringBuilder();
            oechem.OEAddExplicitHydrogens(mol);
            for(StructureAlert alert:alerts){
                oechem.OEPrepareSearch(mol,alert.subSearch);
                int count =0;
                OEMatchBaseIter matches = alert.subSearch.Match(mol);
                for(OEMatchBase ma: matches){
                    if(ma.NumAtoms()==0){
                        continue;
                    }
                    if(!ma.IsValid()){
                        continue;
                    }
                    count += 1;
                }
                if(count <= alert.count){
                    continue;
                }

                boolean matched = false;
                matches.ToFirst();
                for(OEMatchBase ma: matches){
                    if(ma.NumAtoms()==0){
                        continue;
                    }
                    if(!ma.IsValid()){
                        continue;
                    }
                    if (!matched) {
                        matched = true;
                        sb.append(alert);
                        sb.append(";");
                    }
                    for (OEMatchPairAtom atm : ma.GetAtoms()) {
                        if (!structure_alert_set.HasAtom(atm.getTarget())) {
                            structure_alert_set.AddAtom(atm.getTarget());
                        }
                    }
                    for (OEMatchPairBond bnd : ma.GetBonds()) {
                        if (!structure_alert_set.HasBond(bnd.getTarget())) {
                            structure_alert_set.AddBond(bnd.getTarget());
                        }
                    }
                }
            }
            if(sb.length()>0){
                addProperty("Structure Alerts",sb.toString());
                needUpdate = true;
            }
            oechem.OESuppressHydrogens(mol);
        }
    }

    public int getNumAmides(){
        String amideSmarts = "[NX3][CX3](=[OX1])[#6]";
        OESubSearch subSearch = new OESubSearch(amideSmarts);
        oechem.OEPrepareSearch(mol,subSearch);
        int numMatch = 0;
        OEMatchBaseIter matchBases = subSearch.Match(mol, true);
        for(OEMatchBase ma: matchBases){
            numMatch += 1;
        }
        return numMatch;
    }

    public boolean isMarked() {
        return isMarked;
    }

    public void setMarked(boolean marked) {
        isMarked = marked;
    }

    public boolean passBadFragmentsCheck(){
        return ChemFunc.check_library_compounds(this);
    }

    public boolean hasMol3d(){
        return mol3d!=null;
    }


    public OEGraphMol getMol3d(){
        if(mol3d==null){
            mol3d = OEChemFunc.getInstance().getMol3D(mol);
        }
        return mol3d;
    }


    public int getNumAromaticRings(){
        int[] parts = new int[mol.GetMaxAtomIdx()];
        oechem.OEAssignAromaticFlags(mol);
        int numAroRings = oechem.OEDetermineAromaticRingSystems(mol,parts);
        HashMap<Integer,Integer> ringAtmCounts = new HashMap<Integer, Integer>();
        for(int i = 1;i<=numAroRings;i++){
            for(OEAtomBase atm:mol.GetAtoms()){
                if(parts[atm.GetIdx()]==i){
                    if(ringAtmCounts.containsKey(i)){
                        Integer count = ringAtmCounts.get(i) +1;
                        ringAtmCounts.put(i,count);
                    }else{
                        ringAtmCounts.put(i,1);
                    }
                }
            }
        }
        int numOfAroRings = 0;
        for(Integer key:ringAtmCounts.keySet()){
            if(ringAtmCounts.get(key)<=6){
                numOfAroRings += 1;
            }else if(ringAtmCounts.get(key)<=10){
                numOfAroRings += 2;
            }else if(ringAtmCounts.get(key)<=14){
                numOfAroRings += 3;
            }else{
                numOfAroRings += 4;
            }
        }
        return numOfAroRings;

    }

    public double getCLogP(){
        if(!propertyMap.containsKey("CLogP")){
            try {
                ChemFunc.calculateOEProperty(Collections.singletonList(this));
            } catch (Exception e) {
                e.printStackTrace();
                return 0.0;
            }
        }
        return getProperty("CLogP").getValue();
    }

    public double getPSA(){
        if(!propertyMap.containsKey("2d PSA")) {
            ChemFunc.calculateOEPropertySingleMol(this);
        }
        return getProperty("2d PSA").getValue();
    }

    public double getMW(){
        if(!propertyMap.containsKey("molecular weight")){
            ChemFunc.calculateOEPropertySingleMol(this);
        }
        return getProperty("molecular weight").getValue();
    }

    public boolean isEmpty(){
        if(mol==null||mol.NumAtoms()==0){
            return true;
        }
        return false;
    }

    public String getUniqName() {
        return uniq_mol_id;
    }

    public OEAtomBondSet getStructure_alert_set() {
        return structure_alert_set;
    }

    public PropertyMolecule(OEGraphMol mol) {
        this.uniq_mol_id = "pmol_"+mol_id;
        mol_id ++;
        setOEMol(mol);
    }

    public void setOEMol3D(OEGraphMol mol3d){
        if(mol3d!=null) {
            if (mol3d.GetDimension() == 3) {
                this.mol3d = new OEGraphMol(mol3d);
            }
        }
    }

    public void setConf3d(OEGraphMol mol){
        if (mol.GetDimension() == 3) {
            String smiles1 = oechem.OEMolToSmiles(mol);
            String smiles2 = oechem.OEMolToSmiles(this.getMol3d());
            if(smiles1.equals(smiles2)) {
                this.mol3d = new OEGraphMol(mol);
            }
        }
    }

    public void setOEMol(OEGraphMol mol) {
        if(this.mol!=null){
            oechem.OECopySDData(mol,this.mol);
        }
        this.mol = new OEGraphMol(mol);
        if(mol.GetDimension()==3){
            this.mol3d = new OEGraphMol(mol);
            this.mol = OEChemFunc.getInstance().getMol2D(mol);
        }
        this.smiles = ChemFunc.getSmiles(this.mol);

        if(propertyMap==null){
            propertyMap = new LinkedHashMap<String, MolProperty>();
        }else{
            Vector<String> properties = new Vector<>();
            properties.addAll(propertyMap.keySet());
            propertyMap.clear();
            for(String tag:properties){
                if(oechem.OEHasSDData(this.mol, tag)){
                    addProperty(tag,oechem.OEGetSDData(this.mol,tag));
                }
            }
        }
        chemaxon_pKa_dict = new LinkedHashMap<Integer, pKa>();
        moka_pKa_dict = new LinkedHashMap<Integer, pKa>();
        structure_alert_set = new OEAtomBondSet();
        if(this.mol == null){
            this.mol = new OEGraphMol();
            this.name = "None";
        }else{
            if(oechem.OEHasSDData(this.mol,"Name")&&!oechem.OEGetSDData(this.mol, "Name").isEmpty()){
                this.name = oechem.OEGetSDData(this.mol, "Name");
            }else{
                this.name = this.mol.GetTitle();
                if(name==null||name.isEmpty()){
                    name = getUniqName();
                }
                oechem.OEAddSDData(this.mol,"Name",name);
            }
//            long bgn_time = System.nanoTime();
//            svgString = new SVGString(ChemFunc.mol2Svg(this));
            heavyAtomCount = 0;
            for(OEAtomBase atm:mol.GetAtoms()){
                if(atm.GetAtomicNum()>1){
                    heavyAtomCount += 1;
                }
            }
//            long elpased_time = TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - bgn_time);
//            System.out.println(String.format("%s %d %d",mol.GetTitle(),mol.NumAtoms(), elpased_time));
        }
        setNeedUpdate(true);
    }

    public void update2DDepiction(){
        if(mol!=null&&mol.GetDimension()==2){
            svgString = new SVGString(ChemFunc.mol2Svg(this));
        }
    }

    public int getHeavyAtomCount() {
        return heavyAtomCount;
    }

    public void addPropertyObj(String tag, MolProperty property){
        propertyMap.put(tag,property);
    }

    public void addProperty(String tag, String property){
        propertyMap.put(tag,new MolProperty(property));
    }

    public boolean hasProperty(String tag){
        if(propertyMap.containsKey(tag)){
            return true;
        }else{
            return false;
        }
    }

    public MolProperty getProperty(String tag){
        return propertyMap.get(tag);
    }

    public Vector<String> getPropertyNames(){
        Vector<String> properties = new Vector<String>();
        for(String key :propertyMap.keySet()){
            properties.add(key);
        }
        return properties;
    }

    public String getFormula(){
        if(formula.length()==0) {
            oemolostream ofs = new oemolostream();
            ofs.SetFormat(OEFormat.MF);
            ofs.openstring();
            oechem.OEWriteMolecule(ofs, getMol());
            formula = ofs.GetString().split(" ")[0].trim();
        }
        return formula;
    }

    public void addConsensusLogP(){
        if(propertyMap.containsKey("Consensus LogP")){
            return;
        }
        String[] LOGP_NAMES = {"CLogP","XLogP","MoKa LogP","ChemAxon LogP"};
        Vector<Double> values = new Vector<Double>();
        for(String logpName:LOGP_NAMES){
            MolProperty logpProperty = propertyMap.get(logpName);
            if(logpProperty!=null&&logpProperty.isNumerical()){
                values.add(logpProperty.getValue());
            }
        }
        int count = values.size();
        if(count >0){
            double sum = 0.0;
            for(Double v:values){
                sum+=v;
            }
            double avg = sum/count;
            double var = 0.0;
            for(Double v:values){
                var += (v-avg)*(v-avg);
            }
            var = var/count;
            double stddev = Math.sqrt(var);
            MolProperty consensusLogP = new MolProperty(String.format("%5.3f",avg));
            consensusLogP.setStddev(stddev);
            addPropertyObj("Consensus LogP",consensusLogP);
        }
    }

    public OEGraphMol getMol() {
        return mol;
    }

    public String getSmiles(){
        if(smiles==null) {
            smiles = ChemFunc.getSmiles(mol);
        }
        return smiles;
    }

    public SVGString getSvgString() {
        if(needUpdate||svgString==null){
            svgString = new SVGString(ChemFunc.mol2Svg(this));
            needUpdate = false;
        }
        return svgString;
    }

    public void setNeedUpdate(boolean needUpdate) {
        this.needUpdate = needUpdate;
    }

    public String getName() {
        return name;
    }

    public boolean isSelected() {
        return isSelected;
    }

    public void setName(String name) {
        this.name = name;
        mol.SetTitle(name);
        oechem.OESetSDData(mol,"Name",name);
        needUpdate = true;
    }

    public void setIsSelected(boolean isSelected) {
        this.isSelected = isSelected;
    }

    @Override
    public int compareTo(Object o) {
        if(o!=null&&o instanceof PropertyMolecule){
            PropertyMolecule m2 = (PropertyMolecule)o;
            return new Integer(this.getMol().NumAtoms()).compareTo(new Integer(m2.getMol().NumAtoms()));
        }
        return 1;
    }

    @Override
    public boolean equals(Object obj) {
        if(!(obj instanceof PropertyMolecule)){
            return false;
        }else{
            PropertyMolecule p = (PropertyMolecule)obj;
            if(p.getUniqName().equals(getUniqName())){
                return true;
            }
            return false;
        }
    }

    @Override
    public int hashCode() {
        return getUniqName().hashCode();
    }

    public SerializableMol getSerializableMol(){
        OEGraphMol mol = new OEGraphMol(mol3d);
        for(String p:getPropertyNames()){
            oechem.OESetSDData(mol,p,getProperty(p).getProperty());
        }
        String molStr = OEChemFunc.getInstance().getStringFromOEMol(mol, OEFormat.SDF);
        mol.delete();
        return new SerializableMol(molStr,OEFormat.SDF);
    }

    public SerializableMol getSerializableMol2D(){
        OEGraphMol mol = new OEGraphMol(this.mol);
        for(String p:getPropertyNames()){
            oechem.OESetSDData(mol,p,getProperty(p).getProperty());
        }
        String molStr = OEChemFunc.getInstance().getStringFromOEMol(mol, OEFormat.SDF);
        mol.delete();
        return new SerializableMol(molStr,OEFormat.SDF);
    }

    public OEGraphMol getSdfMol(){
        OEGraphMol mol = new OEGraphMol(this.mol);
        for(String p:getPropertyNames()){
            oechem.OESetSDData(mol,p.trim(),getProperty(p).getProperty().trim());
        }
        return mol;

    }

    public String getSdfStr(){
        OEGraphMol mol = new OEGraphMol(this.mol);
        for(String p:getPropertyNames()){
            oechem.OESetSDData(mol,p.trim(),getProperty(p).getProperty().trim());
        }
        String molStr = ChemFunc.getMolString(mol);
        mol.delete();
        return molStr;
    }

    public String getSdfStr3d(){
        OEGraphMol mol = new OEGraphMol(this.getMol3d());
        for(String p:getPropertyNames()){
            oechem.OESetSDData(mol,p,getProperty(p).getProperty());
        }
        String molStr = ChemFunc.getMolString(mol);
        mol.delete();
        return molStr;
    }

    public static void main(String[] args) {
        LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
        Logger rootLogger = loggerContext.getLogger("chemaxon");
        rootLogger.setLevel(Level.OFF);

        Vector<PropertyMolecule> molecules = new Vector<>();
        oemolistream ifs = new oemolistream();
        ifs.open("/home/jfeng/000/20000_enum.sdf");
        boolean isCdx = ifs.GetFormat()==OEFormat.CDX;
        Vector<String> existingTags = new Vector<>();
        OEGraphMol mol = new OEGraphMol();
        int count = 0;
        while (oechem.OEReadMolecule(ifs, mol)) {
            if(isCdx){
                Vector<OEGraphMol> partMols = ChemFunc.splitMol(mol);
                for(OEGraphMol partMol:partMols){
                    molecules.add(new PropertyMolecule(partMol));
                }
                break;
            }
            OESDDataIter oesdDataPairs = oechem.OEGetSDDataPairs(mol);
            while (oesdDataPairs.hasNext()) {
                OESDDataPair p = oesdDataPairs.next();
                String tag = p.GetTag();
                if (!existingTags.contains(tag)) {
                    existingTags.add(tag);
                }
            }
            molecules.add(new PropertyMolecule(mol));
            System.out.println(count ++);
        }
        ifs.close();
    }
}

