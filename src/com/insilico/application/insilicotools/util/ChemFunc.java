package com.insilico.application.insilicotools.util;

import chemaxon.calculations.pka.PKaTrainingResult;
import chemaxon.calculations.pka.PKaTrainingUtil;
import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.jep.function.In;
import chemaxon.marvin.calculations.LogPMethod;
import chemaxon.marvin.calculations.logDPlugin;
import chemaxon.marvin.calculations.logPPlugin;
import chemaxon.marvin.calculations.pKaPlugin;
import chemaxon.marvin.plugin.PluginException;
import chemaxon.struc.*;
import com.insilico.application.insilicotools.InSilicoToolOptions;
import com.insilico.application.insilicotools.data.*;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.alert.AlertRule;
import com.insilico.application.insilicotools.gui.modeling.ModelingException;
import com.google.common.collect.HashBasedTable;
import openeye.oechem.*;
import openeye.oegraphsim.*;
import openeye.oemolprop.*;
import openeye.oemedchem.*;
import org.RDKit.*;
import org.apache.xmlrpc.XmlRpcException;
import org.apache.xmlrpc.client.XmlRpcClient;
import org.apache.xmlrpc.client.XmlRpcClientConfigImpl;

import java.io.*;
import java.net.ConnectException;
import java.net.MalformedURLException;
import java.net.URL;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.apache.commons.codec.binary.Base64;

import javax.vecmath.Matrix3d;
import javastat.util.DataManager;
import org.RDKit.RDKFuncs;
import java.net.HttpURLConnection;
/**
 * Created by jfeng1 on 9/9/15.
 */
// docking_method: {confgen|rigid|mininplace|inplace}

public class ChemFunc {

    private static XmlRpcClient admeClient = null;
    private static XmlRpcClient dockingClient = null;
    public static final String DOCKING_METHOD_FLEXIBLE = "confgen";
    public static final String DOCKING_METHOD_RIGID = "rigid";
    public static final String DOCKING_METHOD_MINIMIZE = "mininplace";
    public static final String DOCKING_METHOD_SCOREONLY = "inplace";
    public static final String DOCKING_MODE_SP = "SP";
    public static final String DOCKING_MODE_XP = "XP";
    public static final HashMap<String,String> PROPERTY_DETAILS;
    public static HashMap<String,String> ADME_URLS;
    public static HashMap<String,String> ADME_KEYPROP;
    static{
        PROPERTY_DETAILS = new HashMap<String, String>();
        PROPERTY_DETAILS.put("ALogD","LogD by Schrodinger");
        PROPERTY_DETAILS.put("CLogP","LogP calculated by RDKit");
        PROPERTY_DETAILS.put("CNS MPO","Pfizer Score (0-6) for CNS delivery, >=4 is desirable.");
        PROPERTY_DETAILS.put("CNS PET MPO","Pfizer Score (0-6) for CNS PET ligand delivery, >=3 is desirable.");
        PROPERTY_DETAILS.put("CNS mTEMPO","Another model for CNS delivery, scaled to 0~6, 6 is the best.");
        PROPERTY_DETAILS.put("No. Aromatic Rings","Number of Aromatic Rings");
        PROPERTY_DETAILS.put("No. Aromatic Ring Systems","Number of Aromatic Ring Systems");
        PROPERTY_DETAILS.put("Ligand Efficiency","measurement of the binding  energy per atom of a ligand to its receptor.");
        PROPERTY_DETAILS.put("Ligand Lipophilic Efficiency","For a given compound LLE is defined as the pIC50 (or pEC50) of interest minus the LogP of the compound");
        PROPERTY_DETAILS.put("Principal Moment of Inertia (NPR1,NPR2)","3D Shape descriptors, a random conformation will be generated automatically if the input molecule does not have 3D coordinates");
        PROPERTY_DETAILS.put("Volume3D","3D Molecular Volume, a random conformation will be generated automatically if the input molecule does not have 3D coordinates");
        PROPERTY_DETAILS.put("SurfaceArea3D","3D Surface area, a random conformation will be generated automatically if the input molecule does not have 3D coordinates");
        PROPERTY_DETAILS.put("PolarSurfaceArea3D","3D polar surface area, a random conformation will be generated automatically if the input molecule does not have 3D coordinates");
        PROPERTY_DETAILS.put("DipoleMoment3D","3D dipole moment, a random conformation will be generated automatically if the input molecule does not have 3D coordinates");
        PROPERTY_DETAILS.put("Ames","A model predicting mutagenic potential of chemical compounds");
        PROPERTY_DETAILS.put("HOMO","Energy of highest occupied molecular orbital");
        PROPERTY_DETAILS.put("LUMO","Energy of lowest unoccupied molecular orbital");
        PROPERTY_DETAILS.put("HBS","Intrinsic hydrogen bond");

        ADME_URLS = new HashMap<>();
//        Plasma Protein Binding https://ppb-service-yatai.cellarity.bentoml.ai
//        Hepatocyte Stability https://hepstability-service-yatai.cellarity.bentoml.ai
//        Blood Stability https://bloodstability-service-yatai.cellarity.bentoml.ai
//        hERG https://herg-service-yatai.cellarity.bentoml.ai
//        Permeability https://perm-public-endpoint-yatai.cellarity.bentoml.ai
//        MoKa https://moka-test2-yatai.cellarity.bentoml.ai

        ADME_URLS.put("Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)","https://perm-public-endpoint-yatai.cellarity.bentoml.ai/predict");
        ADME_URLS.put("Cellarity_hERG_IC50(uM)","https://herg-service-yatai.cellarity.bentoml.ai/predict");
        ADME_URLS.put("hepatocyte stability","https://hepstability-service-yatai.cellarity.bentoml.ai/predict");
        ADME_URLS.put("blood stability","https://bloodstability-service-yatai.cellarity.bentoml.ai/predict");
        ADME_URLS.put("plasma protein unbound(human)%","https://ppb-service-yatai.cellarity.bentoml.ai/predict");

        ADME_KEYPROP = new HashMap<>();
        ADME_KEYPROP.put("Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)","Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)");
        ADME_KEYPROP.put("Cellarity_hERG_IC50(uM)","Cellarity_hERG_IC50(uM)");
        ADME_KEYPROP.put("hepatocyte stability","predicted_label");
        ADME_KEYPROP.put("blood stability","probability_stable");
        ADME_KEYPROP.put("plasma protein unbound(human)%","Cellarity Plasma Human Unbound");

    }

    public static final String[] properties = {
                                                "hydrogen-bond acceptors","hydrogen-bond donors","molecular weight","2d PSA",
                                                "sum of formal charges","number of rings", "No. Aromatic Rings","No. Aromatic Ring Systems","rotatable bonds", "Lipinski violations",
                                                "Ligand Efficiency", "Ligand Lipophilic Efficiency","CLogP","MoKa LogP",
                                                "MoKa LogD",
                                                "MoKa pKa",
                                                "ChemAxon pKa",
                                                "ChemAxon LogD",
                                                "CNS MPO",
                                                "ALogD"

                                                //"HOMO","LUMO","HBS"
                                              };

    public static final String[] adme_models = {

            "Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)",
            "Cellarity_hERG_IC50(uM)",
            "hepatocyte stability",
            "blood stability",
            "plasma protein unbound(human)%"
//                                                "hERG",
//                                                "Papp A->B",
//                                                "General Metabolic Stability: T1/2 (min) (Mouse)"
                                              };

    public static final String[] safety_models = {
     //       "Ames",
            "Structure Alerts"};

    public static final String[] other_models = {
                                                    "Principal Moment of Inertia (NPR1,NPR2)",
//                                                    "Volume3D", "SurfaceArea3D","PolarSurfaceArea3D", "DipoleMoment3D"
                                                };

    public static final String[] TIME_CONSUMING_PROPERTIES_MODELS = {"Principal Moment of Inertia (NPR1,NPR2)","Volume3D", "SurfaceArea3D","PolarSurfaceArea3D", "DipoleMoment3D","HOMO","LUMO","HBS"};

    public final static String[] OE_PROPERTIES = {
            "hydrogen-bond acceptors",
            "hydrogen-bond donors",
            "molecular weight",
            "sum of formal charges",
            "number of rings",
            "rotatable bonds",
            "CLogP",
            "2d PSA"
    };


    public static String getBestLogP(Vector<String> selectedProperties){
        if(selectedProperties==null||selectedProperties.isEmpty()){
            return "CLogP";
        }
        if(selectedProperties.contains("ChemAxon LogP")){
            return "ChemAxon LogP";
        }
        return "CLogP";
    }

    public static final int PORT = 9527;
    public static final int PORT_ADME = 9527;
    public static final int DOCKING_PORT = 9527;
    public static final String ADME_HOST = "10.74.2.128";
    public static final String DOCKING_HOST = "10.74.2.128";
    private static ExecutorService executor = Executors.newFixedThreadPool(10);


    private static XmlRpcClient getADMEClient() throws MalformedURLException {
        if(admeClient==null) {
            admeClient = new XmlRpcClient();
            XmlRpcClientConfigImpl impl = new XmlRpcClientConfigImpl();
            impl.setServerURL(new URL(String.format("http://%s:%d",ADME_HOST,PORT_ADME))); //todo:make sure change to PORT in production version
            admeClient.setConfig(impl);
        }
        return admeClient;
    }

    private static XmlRpcClient getDockingClient() throws MalformedURLException {
        if(dockingClient==null) {
            dockingClient = new XmlRpcClient();
            XmlRpcClientConfigImpl impl = new XmlRpcClientConfigImpl();
            impl.setServerURL(new URL(String.format("http://%s:%d",DOCKING_HOST,DOCKING_PORT)));
//            impl.setServerURL(new URL(String.format("http://%s:%d",DOCKING_HOST,DOCKING_PORT)));
            dockingClient.setConfig(impl);
        }
        return dockingClient;
    }

    public static Object generatePymolSession(String jsonStr) throws MalformedURLException, XmlRpcException {
        Object[] args = new Object[]{jsonStr};
        Object pse = getDockingClient().execute("generatePymolSession", args);
        return pse;
    }

    public static String getCanonicalizedStructure(String molString) throws MalformedURLException, XmlRpcException{
        Object[] args = new Object[]{molString};
        Object new_mol_string = getDockingClient().execute("canonicalizer",args);
        return (String)new_mol_string;
    }

    public static Vector<String> getDockingStructures() throws XmlRpcException, MalformedURLException {
        Vector<String> dockingGrids = new Vector<String>();
        Object[] args = new Object[]{};
        Object resultObj = getDockingClient().execute("getAvailableDockingGrids", args);
        if(resultObj instanceof Object[]){
            Object[] results = (Object[])resultObj;
            for(Object r:results){
                dockingGrids.add((String)r);
            }
        }
        return dockingGrids;
    }

    public static Vector<String> getTemplateStructures() throws XmlRpcException, MalformedURLException {
        Vector<String> templates = new Vector<String>();
        Object[] args = new Object[]{};
        Object resultObj = getDockingClient().execute("getAvailableTemplates", args);
        if(resultObj instanceof Object[]){
            Object[] results = (Object[])resultObj;
            for(Object r:results){
                templates.add((String)r);
            }
        }
        return templates;
    }

    public static String convertNMRLib(String pdbStr) throws MalformedURLException, XmlRpcException {
        Object[] args = new Object[]{pdbStr};
        Object molStrMin = getDockingClient().execute("convert_nmr_lib", args);
        if(molStrMin!=null&&!((String)molStrMin).isEmpty()){
            return (String)molStrMin;
        }
        return null;
    }


    public static String minimize_mol(String molStr) throws MalformedURLException, XmlRpcException {
        Object[] args = new Object[]{molStr};
        Object molStrMin = getDockingClient().execute("minimize_structure", args);
        if(molStrMin!=null&&!((String)molStrMin).isEmpty()){
            return (String)molStrMin;
        }
        return null;
    }

    public static String minimize_mol_constrained(String molStr, String idListStr) throws MalformedURLException, XmlRpcException {
        Object[] args = new Object[]{molStr,idListStr};
        Object molStrMin = getDockingClient().execute("minimize_structure_fixed", args);
        if(molStrMin!=null&&!((String)molStrMin).isEmpty()){
            return (String)molStrMin;
        }
        return null;
    }

/*
    public static String szybki_mol(String molStr) throws MalformedURLException, XmlRpcException {
        Object[] args = new Object[]{molStr};
        Object molStrMin = getDockingClient().execute("szybki", args);
        if(molStrMin!=null&&!((String)molStrMin).isEmpty()){
            return (String)molStrMin;
        }
        return null;
    }

    public static String szybki_mol_constrained(String molStr, String idListStr) throws MalformedURLException, XmlRpcException {
        Object[] args = new Object[]{molStr,idListStr};
        Object molStrMin = getDockingClient().execute("szybki_fixed", args);
        if(molStrMin!=null&&!((String)molStrMin).isEmpty()){
            return (String)molStrMin;
        }
        return null;
    }
*/
    public static String getEvotecInventoryInfo(double cutoff) throws MalformedURLException,XmlRpcException{
        String result = (String)getADMEClient().execute("check_evotec_inventory", new Object[]{new Double(cutoff).toString()});
        return result;
    }

    public static Vector<OEMol> generateConformers(Vector<PropertyMolecule> propertyMolecules, double ewindow, double rms){
        HashMap<String,PropertyMolecule> propertyMolDict = new HashMap<>();
        StringBuilder sb = new StringBuilder();
        for(PropertyMolecule propertyMolecule:propertyMolecules){
            String title = propertyMolecule.getUniqName();
            propertyMolDict.put(title,propertyMolecule);
            OEGraphMol mol = new OEGraphMol(propertyMolecule.getMol3d());
            //Protonator.getInstance().protonate(mol);
            oechem.OEClearSDData(mol);
            mol.SetTitle(title);
            String molStr = OEChemFunc.getInstance().getStringFromOEMol(mol);
            sb.append(molStr);
        }
        Object[] args = new Object[]{
                sb.toString(),
                ewindow,
                rms
        };

        try {
            String resultStr = (String)getDockingClient().execute("oeconfgen",args);
            if(resultStr==null||resultStr.isEmpty()){
                throw new XmlRpcException("Failed to generate conformers...");
            }
            oemolistream ifs = new oemolistream();
            ifs.SetFormat(OEFormat.SDF);
            ifs.openstring(resultStr);
            HashMap<String,OEMol> oemolDict = new HashMap<>();
            Vector<OEMol> oemols = new Vector<>();
            OEGraphMol oeGraphMol = new OEGraphMol();
            while(oechem.OEReadMolecule(ifs,oeGraphMol)){
                String title = oeGraphMol.GetTitle();
                if(!oemolDict.containsKey(title)){
                    OEMol oemol = new OEMol(oeGraphMol);
                    oemolDict.put(title,oemol);
                    PropertyMolecule propertyMol = propertyMolDict.get(title);
                    OEGraphMol tmpMol = propertyMol.getSdfMol();
                    oechem.OECopySDData(oemol,tmpMol);
                    oemols.add(oemol);
                }else{
                    OEMol oemol = oemolDict.get(title);
                    oemol.NewConf(oeGraphMol);
                }
            }
            return oemols;
        } catch (XmlRpcException e) {
            e.printStackTrace();
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }
        return new Vector<>();
    }


    public static Vector<PropertyMolecule> oedock(String receptorName, Vector<PropertyMolecule> inputMols,
                                                String dockingMethod, boolean protonation, String refLigandStr, int numPoses, ProgressReporter progressReporter) throws ModelingException{
        //xmlrpc_docking(self,grid_dir,ligandMolString, poses_per_ligand, precision="SP")
        Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
        Object[] args = new Object[]{receptorName,ChemFunc.convertMolVectorToSDFStringForDocking(inputMols, protonation, progressReporter),
                    dockingMethod, numPoses};
        if(progressReporter!=null) {
            progressReporter.reportProgress("Docking molecules ...", DesignProgressMonitor.INDETERMINATE);
        }
        HashMap<Integer,String> idToNameDict = new HashMap<Integer, String>();
        HashMap<Integer,Vector<PropertyMolecule>> poseDict = new HashMap<Integer, Vector<PropertyMolecule>>();
        int id = 0;
        for(PropertyMolecule m:inputMols){
            String moka_id = oechem.OEGetSDData(m.getMol(),"moka_id");
            idToNameDict.put(id,m.getName());
            m.getMol().SetTitle(moka_id);
            poseDict.put(id,new Vector<PropertyMolecule>());
            id++;
        }
        try {
            String resultStr = (String)getDockingClient().execute("oedocking",args);
            JSONParser parser = new JSONParser();
            JSONObject jsonObject = (JSONObject) parser.parse(resultStr);
            String ligandMolStr = (String)jsonObject.get("ligand");
            if(ligandMolStr.isEmpty()){
                throw new ModelingException("Docking is unable to find a good pose.");
            }
            oemolistream ifs = new oemolistream();
            ifs.SetFormat(OEFormat.SDF);
            ifs.openstring(ligandMolStr);
            OEGraphMol resultMol = new OEGraphMol();
            Vector<Integer> idList = new Vector<Integer>();
            OEGraphMol refMol = null;
            if(refLigandStr!=null){
                refMol = OEChemFunc.getInstance().getOEMolFromString(refLigandStr);
            }
            while(oechem.OEReadMolecule(ifs,resultMol)){
                double comboScoreOE = 0.0;
                if(refMol!=null) {
                    comboScoreOE = OEChemFunc.getInstance().getComboScoreOE(refMol, resultMol);
                }
                PropertyMolecule pmol = new PropertyMolecule(resultMol);
                int moka_id = Integer.parseInt(oechem.OEGetSDData(resultMol,"moka_id"));
                pmol.setName(idToNameDict.get(moka_id));
                poseDict.get(moka_id).add(pmol);
                if(!idList.contains(new Integer(moka_id))){
                    idList.add(moka_id);
                }
                pmol.addProperty("Overlay Score",String.format("%5.2f",comboScoreOE));
                if(dockingMethod.startsWith("posit")){
                    String probability = oechem.OEGetSDData(resultMol, "POSIT::Probability");
                    pmol.addProperty("Probability",probability);
                }else{
                    String chemgaussScore = oechem.OEGetSDData(resultMol, "HYBRID Chemgauss4 score");
                    pmol.addProperty("ChemGauss4 Score",chemgaussScore);
                }
            }
            Collections.sort(idList);
            for(Integer seq:idList){
                Vector<PropertyMolecule> molList = poseDict.get(seq);
                if(!molList.isEmpty()){
                    if(dockingMethod.startsWith("posit")){
                        Collections.sort(molList, new Comparator<PropertyMolecule>() {
                            @Override
                            public int compare(PropertyMolecule o1, PropertyMolecule o2) {
                                double score1 = o1.getProperty("Probability").getValue();
                                double score2 = o2.getProperty("Probability").getValue();
                                return new Double(score2).compareTo(new Double(score1));
                            }
                        });
                    }else{
                        Collections.sort(molList, new Comparator<PropertyMolecule>() {
                            @Override
                            public int compare(PropertyMolecule o1, PropertyMolecule o2) {
                                double score1 = o1.getProperty("ChemGauss4 Score").getValue();
                                double score2 = o2.getProperty("ChemGauss4 Score").getValue();
                                return new Double(score1).compareTo(new Double(score2));
                            }
                        });
                    }
                }
                int n = 0;
                while(n<Math.min(numPoses,molList.size())){
                    propertyMolecules.add(molList.get(n));
                    n++;
                }
            }
            return propertyMolecules;
        } catch (XmlRpcException e) {
            throw(new ModelingException(e.getMessage()));
        } catch (MalformedURLException e) {
            throw(new ModelingException(e.getMessage()));
        } catch (ParseException e) {
            throw(new ModelingException(e.getMessage()));
        }
    }

    public static Vector<PropertyMolecule> rigid_overlay(OEGraphMol reference, Vector<PropertyMolecule> target_mols, ProgressReporter progressReporter){
        Vector<PropertyMolecule> propertyMolecules = new Vector<>();
        Object[] args = new Object[]{OEChemFunc.getInstance().getStringFromOEMol(reference),ChemFunc.convertMolVectorToSDFStringForOverlay(target_mols, false,progressReporter)};
        try {
            if(progressReporter!=null){
                progressReporter.reportProgress("Running Rigid Body Overlaying ...",DesignProgressMonitor.INDETERMINATE);
            }
            String resultStr = (String)getDockingClient().execute("shape_screen_rigid",args);
            oemolistream ifs = new oemolistream();
            ifs.SetFormat(OEFormat.SDF);
            ifs.openstring(resultStr);
            OEGraphMol mol = new OEGraphMol();
            NumberFormat nf = new DecimalFormat("#.##");
            while(oechem.OEReadMolecule(ifs,mol)){
                PropertyMolecule propertyMolecule = new PropertyMolecule(mol);
                int moka_id = Integer.parseInt(oechem.OEGetSDData(mol,"moka_id"));
                double overlayScore = Double.parseDouble(oechem.OEGetSDData(mol,"r_phase_Shape_Sim"));
                propertyMolecule.addProperty("Overlay Score", nf.format(overlayScore));
                propertyMolecules.add(propertyMolecule);
            }
        } catch (XmlRpcException e) {
            e.printStackTrace();
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }

        return propertyMolecules;
    }

    public static Vector<PropertyMolecule> batch_overlay(OEGraphMol reference, Vector<PropertyMolecule> target_mols, double e_window, double rmsd, int nSolutions, ProgressReporter progressReporter){
        Vector<PropertyMolecule> propertyMolecules = new Vector<>();
        Object[] args = new Object[]{OEChemFunc.getInstance().getStringFromOEMol(reference),ChemFunc.convertMolVectorToSDFStringForOverlay(target_mols, false,progressReporter),e_window, rmsd};
        try {
            if(progressReporter!=null){
                progressReporter.reportProgress("Running Overlaying ...",DesignProgressMonitor.INDETERMINATE);
            }
            String resultStr = (String)getDockingClient().execute("shape_screen_flex",args);
            oemolistream ifs = new oemolistream();
            ifs.SetFormat(OEFormat.SDF);
            ifs.openstring(resultStr);
            OEGraphMol mol = new OEGraphMol();
            NumberFormat nf = new DecimalFormat("#.##");
            HashMap<Integer,Vector<PropertyMolecule>> resultMap = new HashMap<>();
            while(oechem.OEReadMolecule(ifs,mol)){
                    PropertyMolecule propertyMolecule = new PropertyMolecule(mol);
                    int moka_id = Integer.parseInt(oechem.OEGetSDData(mol,"moka_id"));
                    double strainEnergy = Double.parseDouble(oechem.OEGetSDData(mol,"r_mmod_Relative_Potential_Energy-S-OPLS"))/4.17;
                    double overlayScore = Double.parseDouble(oechem.OEGetSDData(mol,"r_phase_Shape_Sim"));
                    propertyMolecule.addProperty("Strain Energy", nf.format(strainEnergy));
                    propertyMolecule.addProperty("Overlay Score", nf.format(overlayScore));
                    if(strainEnergy <= e_window) {
                        if (!resultMap.containsKey(new Integer(moka_id))) {
                            resultMap.put(moka_id, new Vector<>());
                        }
                        resultMap.get(new Integer(moka_id)).add(propertyMolecule);
                    }
                    //propertyMolecules.add(propertyMolecule);
            }
            for(Integer key:resultMap.keySet()){
                Collections.sort(resultMap.get(key), new Comparator<PropertyMolecule>() {
                    @Override
                    public int compare(PropertyMolecule o1, PropertyMolecule o2) {
                        double score1 = o1.getProperty("Overlay Score").getValue();
                        double score2 = o2.getProperty("Overlay Score").getValue();
                        return new Double(score2).compareTo(score1);
                    }
                });
                int i = 0;
                for(PropertyMolecule m : resultMap.get(key)){
                    if(i<nSolutions){
                        propertyMolecules.add(m);
                    }
                    i++;
                }
            }
        } catch (XmlRpcException e) {
            e.printStackTrace();
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }
        return propertyMolecules;
    }

//    USE_REF_LIGAND  True
//    REF_LIGAND_FILE refLigandH.sdf
//    CORE_DEFINITION smarts
//    CORE_ATOMS      1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21
//    CORE_SMARTS     O=CNc1cccc(c1)C1CCCN(C1)c1ncnc2[nH]ccc12

//def xmlrpc_docking(self,grid_dir,ligandMolString, poses_per_ligand, precision="SP", use_ref = False, refLigandStr = None,
//                   core_atoms = None, core_smarts = None):
//    result = threads.deferToThread(self.docking, grid_dir,ligandMolString,poses_per_ligand,precision, use_ref,
//    refLigandStr, core_atoms, core_smarts)
//            return result
// docking_method: {confgen|rigid|mininplace|inplace}

    public static Vector<PropertyMolecule> dock(String receptorName, Vector<PropertyMolecule> inputMols, int numPoses,
                                                String dockingMethod, String dockingMode, boolean useRef, String refLigandStr,
                                                String core_atoms, String core_smarts, boolean protonation, ProgressReporter progressReporter) throws ModelingException{
        //xmlrpc_docking(self,grid_dir,ligandMolString, poses_per_ligand, precision="SP")
        Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
        int minPoses = 5;
        Object[] args = null;
        if(useRef&&refLigandStr!=null&&core_atoms!=null&&core_smarts!=null){
            args = new Object[]{receptorName,ChemFunc.convertMolVectorToSDFStringForDocking(inputMols, protonation,progressReporter),
                    numPoses<minPoses?minPoses:numPoses, dockingMethod, dockingMode, true,refLigandStr,core_atoms,core_smarts};
        }else{
            args = new Object[]{receptorName,ChemFunc.convertMolVectorToSDFStringForDocking(inputMols, protonation, progressReporter),
                    numPoses<minPoses?minPoses:numPoses, dockingMethod, dockingMode};
        }
        if(progressReporter!=null) {
            progressReporter.reportProgress("Docking molecules ...", DesignProgressMonitor.INDETERMINATE);
        }
        HashMap<Integer,String> idToNameDict = new HashMap<Integer, String>();
        HashMap<Integer,Vector<PropertyMolecule>> poseDict = new HashMap<Integer, Vector<PropertyMolecule>>();
        int id = 1;
        for(PropertyMolecule m:inputMols){
            idToNameDict.put(id,m.getName());
            poseDict.put(id,new Vector<PropertyMolecule>());
            id++;
        }
        try {
            String resultStr = (String)getDockingClient().execute("docking",args);
            JSONParser parser = new JSONParser();
            JSONObject jsonObject = (JSONObject) parser.parse(resultStr);
            String ligandMolStr = (String)jsonObject.get("ligand");
            if(ligandMolStr.isEmpty()){
                throw new ModelingException("Docking is unable to find a good pose.");
            }
            oemolistream ifs = new oemolistream();
            ifs.SetFormat(OEFormat.SDF);
            ifs.openstring(ligandMolStr);
            OEGraphMol resultMol = new OEGraphMol();
            Vector<Integer> idList = new Vector<Integer>();
            OEGraphMol refMol = null;
            if(refLigandStr!=null){
                refMol = OEChemFunc.getInstance().getOEMolFromString(refLigandStr);
            }
            while(oechem.OEReadMolecule(ifs,resultMol)){
//                double comboScoreOE = 0.0;
//                if(refMol!=null) {
//                     comboScoreOE = OEChemFunc.getInstance().getComboScoreOE(refMol, resultMol);
//                }
                String r_i_docking_score = oechem.OEGetSDData(resultMol, "r_i_docking_score");
                String r_i_glide_einternal = oechem.OEGetSDData(resultMol,"r_i_glide_einternal");
                String i_i_glide_lignum = oechem.OEGetSDData(resultMol,"i_i_glide_lignum");
                String r_i_glide_emodel = oechem.OEGetSDData(resultMol,"r_i_glide_emodel");
                PropertyMolecule pmol = new PropertyMolecule(resultMol);
//                pmol.addProperty("Overlay Score",String.format("%5.2f",comboScoreOE));
                pmol.addProperty("E_Internal",r_i_glide_einternal);
                pmol.addProperty("Interaction Energy",r_i_glide_emodel);
                pmol.addProperty("Docking Score", r_i_docking_score);
                if(r_i_docking_score!=null&&r_i_glide_einternal!=null) {
                    double docking_energy = 0.0;
                    try {
                        docking_energy = Double.parseDouble(r_i_docking_score)/*+Double.parseDouble(r_i_glide_einternal)*/;
                        pmol.addProperty("Glide Score", String.format("%5.2f",docking_energy));
                    } catch (NumberFormatException e) {
                        e.printStackTrace();
                        pmol.addProperty("Glide Score", "0.0");
                    }
                }
                if(i_i_glide_lignum!=null){
                    try {
                        int lignum = Integer.parseInt(i_i_glide_lignum);
                        if(idToNameDict.containsKey(lignum)){
                            pmol.setName(idToNameDict.get(lignum));
                            poseDict.get(lignum).add(pmol);
                            if(!idList.contains(lignum)){
                                idList.add(lignum);
                            }
                        }
                    } catch (NumberFormatException e) {
                        System.err.println(e.getMessage());
                    }
                }
            }
            Collections.sort(idList);
            for(Integer seq:idList){
                Vector<PropertyMolecule> molList = poseDict.get(seq);
                if(!molList.isEmpty()){
                    Collections.sort(molList, new Comparator<PropertyMolecule>() {
                        @Override
                        public int compare(PropertyMolecule o1, PropertyMolecule o2) {
                            double score1 = o1.getProperty("Interaction Energy").getValue();
                            double score2 = o2.getProperty("Interaction Energy").getValue();
                            return new Double(score1).compareTo(new Double(score2));
                        }
                    });
                }
                int n = 0;
                while(n<Math.min(numPoses,molList.size())){
                    propertyMolecules.add(molList.get(n));
                    n++;
                }
            }
            return propertyMolecules;
        } catch (XmlRpcException e) {
            throw(new ModelingException(e.getMessage()));
        } catch (MalformedURLException e) {
            throw(new ModelingException(e.getMessage()));
        } catch (ParseException e) {
            throw(new ModelingException(e.getMessage()));
        }
    }

    public static Vector<String> getReceptorAndLigand(String receptorName) throws XmlRpcException, MalformedURLException, ModelingException {
        Object[] args = new Object[]{receptorName};
        System.out.println(receptorName);
        Object resultObj = getDockingClient().execute("getLigandAndReceptor", args);
        if(resultObj instanceof Object[]){
            Object[] results = (Object[])resultObj;
            if(results.length==2){
                Vector<String> ligandAndReceptor = new Vector<String>();
                ligandAndReceptor.add((String)results[0]);
                ligandAndReceptor.add((String)results[1]);
                return ligandAndReceptor;
            }
        }else{
            if(resultObj instanceof String){
                throw new ModelingException((String)resultObj);
            }
        }
        return null;

    }

    public static String getTemplate(String templateName) throws XmlRpcException, MalformedURLException, ModelingException {
        Object[] args = new Object[]{templateName};
        Object resultObj = getDockingClient().execute("getTemplate", args);
        if(resultObj instanceof String){
            String result = (String)resultObj;
            if(!result.equals("Wrong template name!")){
                return result;
            }else{
                throw new ModelingException(result);
            }
        }else{
                throw new ModelingException("Unknown error.");
        }
    }


    private static double hump(double val, double p0, double p1, double p2, double p3) {
        if (val <= p0 && val > p3) {
            return 0.0;
        } else if (val > p1 && val <= p2) {
            return 1.0;
        } else if (val > p0 && val <= p1) {
            return (val - p0) / (p1 - p0);
        } else if (val > p2 && val <= p3) {
            return (p3 - val) / (p3 - p2);
        } else {
            return 0;
        }
    }

    private static double monotonic(double val, double min, double max) {
        if (val <= min) {
            return 1;
        } else if (val > max) {
            return 0;
        } else {
            return (max - val) / (max - min);
        }
    }

    /*
     1. Select a compound from all example
        compounds and calculate similarity values (Tanimoto coefficient37)
        with the remaining compounds. Then, count the
        number of compounds that have similarity values higher than
        a user-defined similarity cutoff. In this step, two types of
        fingerprints (ECFPs and MDL public keys38) were used for
        structural descriptors to calculate similarity values.

    2. Repeat step 1 for all compounds.

    3. Based on the results of steps 1
       and 2, identify a compound (or compounds) with the largest
       number of similar structures (compounds with Tanimoto
       similarity values higher than the user-defined cutoff value).
       This compound becomes the first key compound candidate.
       If there are multiple compounds with the same number of
       similar structures, then they are all regarded as key compound
       candidates.

    4. Remove the key compound(s) identified in step 3 and their clusters (similar compounds) from the
       population.

    5. Repeat steps 1-4 to find other key compound
       candidates. In this method, steps 1-4 were repeated 5 times
       to identify at least 5 key compound candidates. When the
       predicted key compounds feature several possible stereoisomers,
       they are regarded as one compound.
     */

    public static void predictKeyCompounds(Vector<PropertyMolecule> propertyMolecules, float cutoff, int numHits, int maxRadius,  ProgressReporter progressReporter){
        HashBasedTable<Integer,Integer,Float> tanimotoDict = HashBasedTable.create();
        HashMap<Integer,ExplicitBitVect> fingerPrintHashMap = new HashMap<Integer, ExplicitBitVect>();
        if(progressReporter!=null){
            progressReporter.reportProgress("Generating similairty matrix ...",20);
        }
        for(int i=0;i<propertyMolecules.size();i++){
            PropertyMolecule propertyMolecule = propertyMolecules.get(i);
            org.RDKit.ROMol mol = (org.RDKit.ROMol)org.RDKit.RWMol.MolFromMolBlock(propertyMolecule.getSdfStr());
            ExplicitBitVect fingerprint = RDKFuncs.getMorganFingerprintAsBitVect(mol, 3, 1024);
//            OEFingerPrint fingerPrint = new OEFingerPrint();
//            oegraphsim.OEMakeCircularFP(fingerPrint, propertyMolecule.getMol(), 4096, 0, maxRadius, OEFPAtomType.DefaultCircularAtom, OEFPBondType.DefaultCircularBond);
            fingerPrintHashMap.put(i,fingerprint);
        }
        for(int i=0;i<propertyMolecules.size()-1;i++){
            for(int j=i+1;j<propertyMolecules.size();j++){
                double v = RDKFuncs.TanimotoSimilarity(fingerPrintHashMap.get(i), fingerPrintHashMap.get(j));
                float tanimoto = (float) v;
                tanimotoDict.put(i,j,tanimoto);
                tanimotoDict.put(j,i,tanimoto);
            }
        }
        if(progressReporter!=null){
            progressReporter.reportProgress("Searching key compounds ...",40);
        }
        Set<Integer> bestCompounds = new HashSet<Integer>();
        Set<Integer> neighbors = new HashSet<Integer>();
        for(int k=0;k<numHits;k++) {
            HashMap<Integer, Integer> nbrCountDict = new HashMap<Integer, Integer>();
            HashMap<Integer, Set<Integer>> nbrSetDict = new HashMap<Integer, Set<Integer>>();
            int maxNbrs = -1;
            int keyCmpdId = -1;
            for (int i = 0; i < propertyMolecules.size(); i++) {
                if(neighbors.contains(i)){
                    continue;
                }
                nbrCountDict.put(i, 0);
                nbrSetDict.put(i, new HashSet<Integer>());
                for (int j = 0; j < propertyMolecules.size(); j++) {
                    if(neighbors.contains(j)){
                        continue;
                    }
                    if (j != i) {
                        if (tanimotoDict.get(i, j) > cutoff) {
                            nbrCountDict.put(i, nbrCountDict.get(i) + 1);
                            nbrSetDict.get(i).add(j);
                        }
                    }
                }
                if (nbrCountDict.get(i) > maxNbrs) {
                    maxNbrs = nbrCountDict.get(i);
                    keyCmpdId = i;
                }
            }
            if(keyCmpdId!=-1) {
                System.out.println("KeyCmpdId:"+keyCmpdId);
                System.out.println("No. Nbrs:"+maxNbrs);
                System.out.println("Alt No. Nbrs:"+nbrSetDict.get(keyCmpdId).size());
                bestCompounds.add(keyCmpdId);
                neighbors.addAll(nbrSetDict.get(keyCmpdId));
            }
        }
        if(progressReporter!=null){
            progressReporter.reportProgress("Finishing ...",80);
        }
        for(int i=0;i<propertyMolecules.size();i++){
            if(bestCompounds.contains(i)){
                propertyMolecules.get(i).addProperty("KeyCompound","Yes");
            }else{
                propertyMolecules.get(i).addProperty("KeyCompound","No");
            }
        }
        if(progressReporter!=null){
            progressReporter.reportProgress("Finished",100);
        }
    }

    private static String convertMolVectorToSDFStringForDocking(Vector<PropertyMolecule> mols, boolean protonation, ProgressReporter progressReporter){
        oemolostream ofs = new oemolostream();
        ofs.SetFormat(OEFormat.SDF);
        ofs.openstring();
        int idx = 0;
        int size = mols.size();
        for(PropertyMolecule mol:mols){
            if(progressReporter!=null){
                progressReporter.reportProgress("Generating starting conformation ...",100*idx/size);
            }

            OEGraphMol mol3d = mol.getMol3d();
            if(mol3d==null||mol3d.NumAtoms()==0){
                System.err.println("Failed to generate conformer for "+mol.getSmiles());
                idx++;
                continue;
            }
            OEGraphMol oemol = new OEGraphMol(mol3d);
            oemol.SetTitle(mol.getName());
            OEGroupBaseIter iter = oemol.GetGroups(new OEIsMDLStereoGroup());
            for(OEGroupBase g :iter) {
                oemol.DeleteGroup(g);
            }
            if (protonation) {
                Protonator.getInstance().protonate(oemol);
            } else {
                Protonator.getInstance().adjust_hydrogens(oemol, false);
            }

            String molStr = null;
            try {
                molStr = minimize_mol(OEChemFunc.getInstance().getStringFromOEMol(oemol));
                oemol = OEChemFunc.getInstance().getOEMolFromString(molStr);
            } catch (MalformedURLException e) {
                e.printStackTrace();
            } catch (XmlRpcException e) {
                e.printStackTrace();
            }
            //oemol = OEChemFunc.getInstance().minimizeOEMol(oemol,null,protonation);
            oechem.OEClearSDData(oemol);
            oechem.OESetSDData(oemol, "moka_id", String.format("%d",idx++));
            oechem.OEWriteMolecule(ofs, oemol);
            oemol.delete();
        }
        return ofs.GetString();
    }

    private static String convertMolVectorToSDFStringForOverlay(Vector<PropertyMolecule> mols, boolean protonation, ProgressReporter progressReporter){
        oemolostream ofs = new oemolostream();
        ofs.SetFormat(OEFormat.SDF);
        ofs.openstring();
        int idx = 0;
        int size = mols.size();
        for(PropertyMolecule mol:mols){
            if(progressReporter!=null){
                progressReporter.reportProgress("Generating starting conformation ...",100*idx/size);
            }

            OEGraphMol mol3d = mol.getMol3d();
            if(mol3d==null||mol3d.NumAtoms()==0){
                System.err.println("Failed to generate conformer for "+mol.getSmiles());
                idx++;
                continue;
            }
            OEGraphMol oemol = new OEGraphMol(mol3d);
            if (protonation) {
                Protonator.getInstance().protonate(oemol);
            } else {
                Protonator.getInstance().adjust_hydrogens(oemol, false);
            }

            //oemol = OEChemFunc.getInstance().minimizeOEMol(oemol,null,protonation);
            oechem.OESetSDData(oemol, "moka_id", String.format("%d",idx++));
            oechem.OEWriteMolecule(ofs, oemol);
            oemol.delete();
        }
        return ofs.GetString();
    }


    private static String convertMolVectorToSDFString(Vector<PropertyMolecule> mols){
        oemolostream ofs = new oemolostream();
        ofs.SetFormat(OEFormat.SDF);
        ofs.openstring();
        int idx = 0;
        for(PropertyMolecule mol:mols){
            OEGraphMol oemol = new OEGraphMol(mol.getMol());
            OEGroupBaseIter iter = oemol.GetGroups(new OEIsMDLStereoGroup());
            for(OEGroupBase g :iter) {
                oemol.DeleteGroup(g);
            }
            oechem.OEClearSDData(oemol);
            for(String p:ChemFunc.OE_PROPERTIES){
                if(mol.getPropertyNames().contains(p)) {
                    oechem.OESetSDData(oemol, p, mol.getProperty(p).getProperty());
                }
            }
            oemol.SetTitle(String.format("%d",idx));
            oechem.OESetSDData(oemol, "moka_id", String.format("%d",idx++));
            oechem.OEWriteMolecule(ofs, oemol);
            oemol.delete();
        }
        return ofs.GetString();
    }

    private static String convertMolVectorToSDFString3D(Vector<PropertyMolecule> mols){
        oemolostream ofs = new oemolostream();
        ofs.SetFormat(OEFormat.SDF);
        ofs.openstring();
        int idx = 0;
        for(PropertyMolecule mol:mols){
            OEGraphMol oemol = new OEGraphMol(mol.getMol3d());
            oechem.OEClearSDData(oemol);
            oechem.OESetSDData(oemol, "moka_id", String.format("%d",idx++));
            oechem.OEWriteMolecule(ofs, oemol);
            oemol.delete();
        }
        return ofs.GetString();
    }

    private static String convertMolVectorToSDFStringEfflux(Vector<PropertyMolecule> mols) throws MalformedURLException, ParseException, ConnectException, XmlRpcException {
        oemolostream ofs = new oemolostream();
        ofs.SetFormat(OEFormat.SDF);
        ofs.openstring();
        int idx = 0;
        ChemFunc.calculateChemAxon_pKa(mols,null);
        ChemFunc.calculateChemAxonLogD(mols,7.4,null);

        for(PropertyMolecule mol:mols){
            OEGraphMol oemol = new OEGraphMol(mol.getMol());
            OEGroupBaseIter iter = oemol.GetGroups(new OEIsMDLStereoGroup());
            for(OEGroupBase g :iter) {
                oemol.DeleteGroup(g);
            }
            oechem.OEClearSDData(oemol);
            for(String p:ChemFunc.OE_PROPERTIES){
                if(mol.getPropertyNames().contains(p)) {
                    oechem.OESetSDData(oemol, p, mol.getProperty(p).getProperty());
                }
            }
            oechem.OESetSDData(oemol,"CLogP",mol.getProperty("ChemAxon LogP")+"");
            MolProperty chemAxon_logD = mol.getProperty("ChemAxon LogD");
            oechem.OESetSDData(oemol,"ChemAxon LogD", chemAxon_logD==null?"0":chemAxon_logD.getProperty());
            MolProperty chemAxon_acidic_pKa = mol.getProperty("ChemAxon Acidic pKa");
            oechem.OESetSDData(oemol,"ChemAxon Acidic pKa", chemAxon_acidic_pKa==null?"0":chemAxon_acidic_pKa.getProperty());
            MolProperty chemAxon_basic_pKa = mol.getProperty("ChemAxon Basic pKa");
            oechem.OESetSDData(oemol,"ChemAxon Basic pKa", chemAxon_basic_pKa==null?"0":chemAxon_basic_pKa.getProperty());
            oechem.OESetSDData(oemol, "moka_id", String.format("%d",idx++));
            oechem.OEWriteMolecule(ofs, oemol);
            oemol.delete();
        }
        return ofs.GetString();
    }


    public static void calculatePMIBatch(Vector<PropertyMolecule> molecules){
        if(molecules!=null&&molecules.size()>0){
            for(PropertyMolecule pm:molecules){
                OEGraphMol mol3d = pm.getMol3d();
                if(mol3d!=null) {
                    double[] pmi = calculatePMI(mol3d);
                    pm.addProperty("NPR1", String.format("%5.2f", pmi[0]));
                    pm.addProperty("NPR2", String.format("%5.2f", pmi[1]));
                }else{
                    pm.addProperty("NPR1", "0.0");
                    pm.addProperty("NPR2", "0.0");
                }
            }
        }
    }

    public static double[] calculatePMI(OEGraphMol mol){
        double[] pmi = new double[2];
        pmi[0] = 0.0;
        pmi[1] = 0.0;
        OEGraphMol mol3d = OEChemFunc.getInstance().getMol3D(mol);
        if(mol3d==null){
            return pmi;
        }
        DataManager dataManager = new DataManager();
        if(mol3d!=null){
            OEFloatArray centerOfMass = new OEFloatArray(new float[3]);
            oechem.OEGetCenterOfMass(mol3d,centerOfMass,true);
            double cmx = centerOfMass.getItem(0);
            double cmy = centerOfMass.getItem(1);
            double cmz = centerOfMass.getItem(2);
            Matrix3d inertia = new Matrix3d(0,0,0,0,0,0,0,0,0);
            for(OEAtomBase atm:mol3d.GetAtoms()){
                double[] xyz = new double[3];
                mol3d.GetCoords(atm,xyz);
                double x = xyz[0]-cmx;
                double y = xyz[1]-cmy;
                double z = xyz[2]-cmz;
                int element_no = atm.GetAtomicNum();
                double atomic_weight = oechem.OEGetAverageWeight(element_no);
                Matrix3d m = new Matrix3d(y*y+z*z,-x*y,-x*z,-y*x,x*x+z*z,-y*z,-z*x,-z*y,x*x+y*y);
                m.mul(atomic_weight);
                inertia.add(m);
            }
            double[][] inertiaMatrix = new double[3][3];
            for(int i=0;i<3;i++){
                for(int j=0;j<3;j++){
                    inertiaMatrix[i][j] = inertia.getElement(i,j);
                }
            }
            double[] eigenvalues = dataManager.eigenvalue(inertiaMatrix);
            Arrays.sort(eigenvalues);
            if(eigenvalues.length==3&&eigenvalues[2]>0){
                pmi[0] = eigenvalues[0]/eigenvalues[2];
                pmi[1] = eigenvalues[1]/eigenvalues[2];
            }
        }
        return pmi;
    }

    public static String formatCYNumber(String rawCYNumber){
        if(rawCYNumber==null||rawCYNumber.isEmpty()){
            return null;
        }
        String rawCYNumber2 = rawCYNumber.trim().toUpperCase();
        if(rawCYNumber2.length()>10){
            return null;
        }
        if(rawCYNumber2.matches("CY-\\d+")){
            if(rawCYNumber2.length()==10) {
                return rawCYNumber2;
            }else{
                String tmp = rawCYNumber2.replaceAll("CY-","");
                if(tmp.matches("\\d+")&&tmp.length()<=6){
                    String rawNumber3 = String.format("CY-2%06d",Integer.parseInt(tmp));
                    return rawNumber3;
                }
            }
        }
        if(rawCYNumber2.matches("\\d+")&&rawCYNumber2.length()<=6){
            String rawNumber3 = String.format("CY-2%06d",Integer.parseInt(rawCYNumber2));
            return rawNumber3;
        }
        return null;
    }

    static final HashMap<String, double[]> ruleDict = new HashMap<String, double[]>() {
        {
            put("ChemAxon LogP", new double[]{3,5});
            put("2d PSA", new double[]{20, 40, 90, 120});
            put("ChemAxon Basic pKa",new double[]{8,10});
            put("hydrogen-bond donors", new double[]{0.5, 3.5});
            put("ChemAxon LogD", new double[]{2,4});
            put("molecular weight", new double[]{360, 500});
        }
    };

    static final HashMap<String, double[]> pet_ruleDict = new HashMap<String, double[]>() {
        {
            put("ChemAxon LogP", new double[]{2.8,4});
            put("2d PSA", new double[]{32.3, 44.8, 63.3, 86.2});
            put("ChemAxon Basic pKa",new double[]{7.2,9.5});
            put("hydrogen-bond donors", new double[]{1, 2});
            put("ChemAxon LogD", new double[]{1.7,2.8});
            put("molecular weight", new double[]{305.3, 350.5});
        }
    };


    public static void generateCNSMPODescriptors(Vector<PropertyMolecule> mols, boolean pet){
        for(PropertyMolecule mol:mols){
            String logP_property;
            String logD_property;
            String basicPka_property;
            MolProperty logp;
            MolProperty logd74;
            MolProperty basicPka;
            logP_property = "ChemAxon LogP";
            logD_property = "ChemAxon LogD";
            basicPka_property = "ChemAxon Basic pKa";
            logp = mol.getProperty(logP_property);
            logd74 = mol.getProperty(logD_property);
            basicPka = mol.getProperty(basicPka_property);

            MolProperty hbd = mol.getProperty("hydrogen-bond donors");
            MolProperty psa = mol.getProperty("2d PSA");
            MolProperty mw = mol.getProperty("molecular weight");
            if(logp!=null&&logp.isNumerical()&&logd74!=null&&logd74.isNumerical()&&hbd!=null&&hbd.isNumerical()&&psa!=null&&psa.isNumerical()&&mw!=null&&mw.isNumerical()){
                HashMap<String,Double> valueDict = new HashMap<String, Double>();
                double basic_pka = 0.0;
                if(basicPka!=null&&basicPka.isNumerical()){
                    basic_pka = basicPka.getValue();
                }
                valueDict.put(logP_property, logp.getValue());
                valueDict.put(logD_property, logd74.getValue());
                valueDict.put(basicPka_property,basic_pka);
                valueDict.put("hydrogen-bond donors",hbd.getValue());
                valueDict.put("2d PSA", psa.getValue());
                valueDict.put("molecular weight",mw.getValue());
                double sum = 0.0;
                String[] propertySet = {logP_property, "2d PSA", basicPka_property, "hydrogen-bond donors", logD_property, "molecular weight"};
                for(String property : propertySet){
                    double value = valueDict.get(property);
                    double[] rules = null;
                    if(pet){
                        rules = pet_ruleDict.get(property);
                    }else{
                        rules = ruleDict.get(property);
                    }
                    if (rules.length == 2) {
                        sum += monotonic(value, rules[0], rules[1]);
                    } else {
                        sum += hump(value, rules[0], rules[1], rules[2], rules[3]);
                    }
                }
                if(pet){
                    mol.addPropertyObj("CNS PET MPO", new CNSMPObj(String.format("%5.2f", sum)));
                }else {
                    mol.addPropertyObj("CNS MPO", new CNSMPObj(String.format("%5.2f", sum)));
                }
            }
        }
    }

    public static void generateNoAromaticRingsDescriptors(Vector<PropertyMolecule> mols){
        for(PropertyMolecule mol:mols){
            OEGraphMol oemol = mol.getMol();
            int[] parts = new int[oemol.GetMaxAtomIdx()];
            oechem.OEFindRingAtomsAndBonds(oemol);
            oechem.OEAssignAromaticFlags(oemol);
            int n = oechem.OEDetermineAromaticRingSystems(oemol, parts);
            mol.addProperty("No. Aromatic Ring Systems", ""+n);
        }
    }

    public static void calculateNoAroRings(Vector<PropertyMolecule> mols){
        for(PropertyMolecule mol:mols){
            ROMol romol = (org.RDKit.ROMol)org.RDKit.RWMol.MolFromMolBlock(mol.getSdfStr());
            int numAromaticRings = (int)RDKFuncs.calcNumAromaticRings(romol);
            mol.addProperty("No. Aromatic Rings",""+numAromaticRings);
        }
    }
    /*
    public static void generateCLogPSingleMol(PropertyMolecule mol) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        Object[] args = new Object[]{mol.getSmiles()};
        String clogp = (String)admeClient.execute("CLogP", args);
        System.out.println(String.format("%s %s",mol.getSmiles(),clogp));
        mol.addProperty("CLogP",clogp);
    }

    public static void generateCLogP(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        StringBuilder sb = new StringBuilder();
        for(PropertyMolecule mol:mols){
            sb.append(String.format("%s %s\n",mol.getSmiles(),mol.getUniqName()));
        }
        Object[] args = new Object[]{sb.toString()};
        String clogpResult = (String)admeClient.execute("CLogPBatch", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(clogpResult);
        for(PropertyMolecule mol:mols){
            if(obj.containsKey(mol.getUniqName())){
                mol.addProperty("CLogP",(String)obj.get(mol.getUniqName()));
            }
        }
    }
     */

    public static void calculateFreeForm(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        dockingClient = getDockingClient();
        String inputSdf = convertMolVectorToSDFString3D(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String freeformResult = (String)dockingClient.execute("freeform", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(freeformResult);
        int idx = 0;
        String[] propertyKeys = new String[]{"Erel","deltaG","Local Strain","Global Strain"};
        for(PropertyMolecule mol:mols){
            //{"ElocStrainFreeSheff":"1.7858","GlblStrainSheff":"4.80","conf_dGSheff":"3.0100"}
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                JSONObject obj1 = (JSONObject)obj.get(moka_id);
                for(String property:propertyKeys) {
                    if (obj1.containsKey(property)) {
                        mol.addProperty(property,(String)obj1.get(property));
                    }
                }
            }
        }
    }

    public static void generateVDss(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();

        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String vdssResult = (String)admeClient.execute("predictVdss", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(vdssResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                mol.addProperty("Human VDss(L/kg)",(obj.get(moka_id)).toString());
            }
        }
    }

    public static void generateInSilicoADMEproperties(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        String [] properties = new String[]{"cHLM","cHPPB","cMDR1","cRLM","cRPPB","cSolubility"};
        admeClient = getADMEClient();

        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String admeResult = (String)admeClient.execute("inSilicoADME", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(admeResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                JSONObject dict = (JSONObject)obj.get(moka_id);
                for(String property:properties) {
                    mol.addProperty(property, ((Double)dict.get(property)).toString());
                }
            }
        }
    }

    public static void generateBBB(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();

        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String vdssResult = (String)admeClient.execute("predictBBB", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(vdssResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                String property = (obj.get(moka_id)).toString();
                double v = Double.parseDouble(property);
                double bbb = Math.pow(10,v);
                mol.addProperty("BB Ratio", nf.format(bbb));
            }
        }
    }

    public static void calculateSolubility_pH_7(Vector<PropertyMolecule> mols)throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        String property = "VolSurf_LgS7";
        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String volsurfResult = (String)admeClient.execute("VolSurfDescriptor", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(volsurfResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                JSONObject dict = (JSONObject) obj.get(moka_id);
                if(dict.containsKey(property)){
                    Object value = dict.get(property);
                    if(value instanceof Double){
                        double new_value = (Double)value;
                        MolProperty p = new MolProperty(nf.format(new_value));
                        mol.addPropertyObj("Log_Solubility_pH7",p);
                    }
                }
            }
        }
    }

    public static void calculateMetabolicStability(Vector<PropertyMolecule> mols)throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        String property = "VolSurf_MetStab";
        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String volsurfResult = (String)admeClient.execute("VolSurfDescriptor", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(volsurfResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                JSONObject dict = (JSONObject) obj.get(moka_id);
                if(dict.containsKey(property)){
                    Object value = dict.get(property);
                    if(value instanceof Double){
                        MolProperty p = new MolProperty(nf.format((Double)value));
                        mol.addPropertyObj("Metabolic Stability",p);
                    }
                }
            }
        }
    }

    public static void generateVolSurfDescriptors(Vector<PropertyMolecule> mols, Vector<String> selectedProperties) throws MalformedURLException, XmlRpcException, ParseException, ConnectException {
        admeClient = getADMEClient();
        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String volsurfResult = (String)admeClient.execute("VolSurfDescriptor", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(volsurfResult);
        int idx = 0;
        if(selectedProperties==null){
            selectedProperties = new Vector<String>();
        }
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                JSONObject dict = (JSONObject) obj.get(moka_id);
                for(Object key : dict.keySet()){
                    Object value = dict.get(key);
                    if(idx==1){
                        if(!selectedProperties.contains((String)key)){
                            selectedProperties.add((String)key);
                        }
                    }
                    if(value instanceof Double){
                        mol.addProperty((String)key,nf.format(value));
                    }else if(value instanceof String) {
                        mol.addProperty((String) key, (String) dict.get(key));
                    }else {
                        mol.addProperty((String)key,(dict.get(key)).toString());
                    }
                }
            }
        }
    }

    public static String getMolString(OEGraphMol mol){
        if(mol!=null){
            oemolostream ofs = new oemolostream();
            ofs.SetFormat(OEFormat.SDF);
            ofs.openstring();
            oechem.OEWriteMolecule(ofs,mol);
            return ofs.GetString();
        }
        return null;
    }

    public static void predictClearanceMechanism(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();

        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String clearanceResult = (String)admeClient.execute("predictClearanceMode", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(clearanceResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                Object o = obj.get(moka_id);
                JSONObject childObj = (JSONObject)parser.parse(o.toString());
                mol.addProperty("Clearance Route", (String)childObj.get("route"));
                mol.addProperty("PC_1",childObj.get("pc_1").toString());
                mol.addProperty("PC_2",childObj.get("pc_2").toString());
            }
        }
    }
//    def xmlrpc_predictMetabolicClearance(self,inputSdf):
//            return self.predictClearance(inputSdf, True)
//
//    def xmlrpc_predictRenalClearance(self,inputSdf):
//            return self.predictClearance(inputSdf,False)

    public static void getMopacProperties(Vector<PropertyMolecule> mols, ProgressReporter progressReporter) throws MalformedURLException, XmlRpcException, ParseException {
        if(mols!=null&&!mols.isEmpty()){
            int numMols = mols.size();
            int count = 0;
            for(PropertyMolecule mol:mols){
                if(progressReporter!=null){
                    progressReporter.reportProgress("Progress",100*++count/numMols);
                }
                OEGraphMol mol3d = mol.getMol3d();
                oechem.OEAssignFormalCharges(mol3d);
                int formal_charges = 0;
                for(OEAtomBase atom: mol3d.GetAtoms()){
                    formal_charges += atom.GetFormalCharge();
                }
                String result = (String)getDockingClient().execute("mopac", new Object[]{mol.getSdfStr3d(), formal_charges});
                JSONParser parser = new JSONParser();
                JSONObject obj = (JSONObject)parser.parse(result);
                float homo = Float.parseFloat((String)obj.get("HOMO"));
                float lumo = Float.parseFloat((String)obj.get("LUMO"));
                System.out.println(homo+" "+lumo);
                JSONObject dict = (JSONObject) obj.get("SAR");

                double o_sum_donors = 0.0;
                double o_sum_acceptors = 0.0;
                double n_sum_donors = 0.0;
                double n_sum_acceptors = 0.0;

                for(Object key:dict.keySet()){
                    JSONObject atomDict = (JSONObject)dict.get(key);
                    int atom_idx = Integer.parseInt((String)key)-1;
                    OEAtomBase atm = mol3d.GetAtom(new OEHasAtomIdx(atom_idx));
                    boolean is_donor = false;
                    for(OEAtomBase nbr:atm.GetAtoms()){
                        if(nbr.IsHydrogen()){
                            is_donor = true;
                            break;
                        }
                    }
                    double se = Double.parseDouble((String)atomDict.get("SE"));
                    double sn = Double.parseDouble((String)atomDict.get("SN"));
                    double asp = Double.parseDouble((String)atomDict.get("ASP"));
                    if(atm.IsOxygen()&&is_donor){
                        o_sum_donors += se;
                    }
//                    if(atm.IsOxygen()){
//                        o_sum_acceptors+=sn;
//                    }
                    if(atm.IsNitrogen()&&is_donor){
                        n_sum_donors += se;
                    }
//                    if(atm.IsNitrogen()&&!is_donor){
//                        n_sum_acceptors += sn;
//                    }
                }
                mol.addProperty("HBS",""+(o_sum_donors+o_sum_acceptors+n_sum_acceptors+n_sum_donors));
                mol.addProperty("HOMO",""+homo);
                mol.addProperty("LUMO",""+lumo);
            }
        }
    }

    public static void predictClearance(Vector<PropertyMolecule> mols, boolean isMetabolic) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        //"Metabolic Clearance","Renal Clearance"
        String functionName = "predictRenalClearance";
        String tagName = "Human Renal Clearance(mL/min/kg)";
        if(isMetabolic){
            functionName = "predictMetabolicClearance";
            tagName = "Human Metabolic Clearance(mL/min/kg)";
        }

        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String clearanceResult = (String)admeClient.execute(functionName, args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(clearanceResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                mol.addProperty(tagName,(obj.get(moka_id)).toString());
            }
        }
    }
    public static void predict_hERG_DL(Vector<PropertyMolecule> mols)throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        //String property = "hERG IC50 Pharmaron_prediction";
        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String result = (String)admeClient.execute("hERG", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(result);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                if(obj.containsKey(moka_id)){
                    double value = (Double) obj.get(moka_id);
                    if(value<=0){
                        value = 0.001;
                    }
                    MolProperty p = new MolProperty(nf.format(value));
                    mol.addPropertyObj("hERG",p);
                }
            }
        }
    }
    public static void predict_MDCK_DL(Vector<PropertyMolecule> mols)throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        String property = "Papp A-B_prediction";
        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String result = (String)admeClient.execute("mdck", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(result);
        System.out.println(result);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                double value = (Double) obj.get(moka_id);
                if(value<=0){
                    value = 0.001;
                }
                MolProperty p = new MolProperty(nf.format(value));
                mol.addPropertyObj("Papp A->B",p);
            }
        }
    }

    public static void predict_MSM_DL(Vector<PropertyMolecule> mols)throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        NumberFormat nf = new DecimalFormat("#.###");
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String result = (String)admeClient.execute("msm", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(result);
        System.out.println(result);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                double value = (Double) obj.get(moka_id);
                if(value<=0){
                    value = 0.001;
                }
                MolProperty p = new MolProperty(nf.format(value));
                mol.addPropertyObj("General Metabolic Stability: T1/2 (min) (Mouse)",p);
            }
        }
    }

    //todo:
    //todo:New ADME models by ML group@Cellarity
    //todo:
    //
    //
//    def calculate_hERG(self, smilesDict):
//    conn = httplib.HTTPConnection("herg.cellarity.int", port=80)
//    headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
//    parser = JSONParser()
//
//    smilesList = JSONArray()
//        for key in smilesDict.keys():
//    smilesEntry = JSONObject()
//            smilesEntry.put("id", key)
//            smilesEntry.put("smiles", smilesDict[key])
//            smilesList.add(smilesEntry)
//            conn.request("POST", "/predict", smilesList.toJSONString(), headers)
//    response = conn.getresponse()
//            if response.status == 200:
//    data = response.read()
//    result_dict = parser.parse(data)
//            return result_dict
    public static void runMLModels(Vector<PropertyMolecule> propertyMolecules, String modelName) throws IOException, ParseException {
        String url = ADME_URLS.get(modelName);
        URL obj = new URL(url);
        HttpURLConnection con = (HttpURLConnection) obj.openConnection();
        con.setRequestProperty("accept","application/json");
        con.setRequestProperty("Content-Type", "application/json");
        con.setRequestMethod("POST");

        JSONParser parser = new JSONParser();
        JSONArray smilesList = new JSONArray();
        int index = -1;
        for(PropertyMolecule propertyMolecule:propertyMolecules) {
            index += 1;
            JSONObject smilesObj = new JSONObject();
            smilesObj.put("id", index);
            smilesObj.put("smiles", oechem.OEMolToSmiles(propertyMolecule.getMol()));
            smilesList.add(smilesObj);
        }
        String parmas = smilesList.toJSONString();

        // For POST only - START
        con.setDoOutput(true);
        OutputStream os = con.getOutputStream();
        os.write(parmas.getBytes());
        os.flush();
        os.close();
        // For POST only - END

        int responseCode = con.getResponseCode();
        //System.out.println("POST Response Code :: " + responseCode);

        if (responseCode == HttpURLConnection.HTTP_OK) { //success
            BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
            String inputLine;
            StringBuffer response = new StringBuffer();

            while ((inputLine = in.readLine()) != null) {
                response.append(inputLine);
            }
            in.close();

            // print result
            JSONObject result = (JSONObject) parser.parse(response.toString());
            for(Object key : result.keySet()){
                int compound_id = Integer.parseInt((String) key);
                JSONObject molResult = (JSONObject) result.get(key);
                for(Object property: molResult.keySet()){
                    String propertyName = (String) property;
                    Object propertyValueObj = molResult.get(property);
                    String propertyValue = null;
                    System.out.println(compound_id+" "+property+"  "+propertyValueObj);
                    if(propertyValueObj instanceof String) {
                        propertyValue = (String) propertyValueObj;
                    }
                    if(propertyValueObj instanceof Long) {
                        propertyValue = Long.toString((Long)propertyValueObj);
                    }

                    if(propertyValueObj instanceof Integer){
                        propertyValue = Integer.toString((Integer)propertyValueObj);
                    }

                    if(propertyValueObj instanceof Float){
                        propertyValue = Float.toString((Float)propertyValueObj);
                    }

                    if(propertyValueObj instanceof Double){
                        propertyValue = Double.toString((Double)propertyValueObj);
                    }
                    if(propertyValue!=null && propertyName.equals(ADME_KEYPROP.get(modelName))) {
                        propertyMolecules.get(compound_id).addProperty(modelName, propertyValue);
                    }else{
                        System.out.println(property);
                    }
                }
            }

        } else {
            throw new IOException("Failed to calculate ");
        }
    }

    public static void generateAlogD(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException {
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFString(mols);
        Object[] args = new Object[]{inputSdf};
        String result = (String)admeClient.execute("alogd",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(result);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)) {
                double value = (Double)obj.get(moka_id);
                mol.addProperty("ALogD", ""+value);
            }
        }
    }

    public static void generateMoKaDescriptors(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException {
        admeClient = getADMEClient();

        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String mokaResult = (String)admeClient.execute("Moka", args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(mokaResult);
        int idx = 0;
        for(PropertyMolecule mol:mols){
            String moka_id = ""+idx++;
            if(obj.containsKey(moka_id)){
                JSONObject o1 = (JSONObject)obj.get(moka_id);
                if(o1.containsKey("MoKa_LogD7.4")){
                    mol.addProperty("MoKa LogD",(String)o1.get("MoKa_LogD7.4"));
                }
                if(o1.containsKey("MoKa_LogP")){
                    mol.addProperty("MoKa LogP",((String)o1.get("MoKa_LogP")));
                }
                if(o1.containsKey("pKa")){
                    JSONArray pkas = (JSONArray) o1.get("pKa");
                    Vector<pKa> aPkas = new Vector<pKa>();
                    Vector<pKa> bPkas = new Vector<pKa>();
                    for(int i=0;i<pkas.size();i++){
                        JSONObject o = (JSONObject) pkas.get(i);
                        float value = ((Double) o.get("value")).floatValue();
                        int atom_id = ((Long) o.get("atomid")).intValue()-1;
                        String t = (String) o.get("type");
                        pKa_type type = pKa_type.BASIC;
                        if(t.equals("a")){
                            type = pKa_type.ACIDIC;
                        }

                        pKa p = new pKa(atom_id,mol,value,type);
                        if(p.isAcidic()){
                            aPkas.add(p);
                        }else{
                            bPkas.add(p);
                        }
                    }
                    if(aPkas.size()>0) {
                        Collections.sort(aPkas);
                        pKa aPKa = aPkas.get(0);
                        float value = aPKa.getValue();
                        mol.addProperty("MoKa Acidic pKa",""+ value);
                        for(pKa p:aPkas) {
                            mol.addPKa(p.getAtom_id(), p, PropertyMolecule.MOKA_PKA);
                        }
                    }
                    if(bPkas.size()>0){
                        Collections.sort(bPkas);
                        float value = bPkas.get(0).getValue();
                        mol.addProperty("MoKa Basic pKa",""+ value);
                        for(pKa p:bPkas){
                            mol.addPKa(p.getAtom_id(),p, PropertyMolecule.MOKA_PKA);
                        }
                    }
                }
            }
        }
//        print s.Moka(ligand)
//        print s.hERG(ligand)
//        print s.PPB(ligand)
    }

//    public static final String[] ADME_MODELS = {"Passive Permeability","PgP Substrate", "Plasma Protein Binding (Human)",
//            "pH-Solubility Profile","pH-Lipophilicity Profile","VolSuf ADME Report","Structure Alerts","hERG inhibition"};

    public static void calculateChemAxon_pKa(Vector<PropertyMolecule> mols, ProgressReporter progressReporter){
        pKaPlugin plugin = new pKaPlugin();

        plugin.setDoublePrecision(2);
        if(InSilicoToolOptions.use_trained_chemaxon_pKa) {
            try {
                URL url = new URL("http://javelin.corp.My.com:8080/insilico/chemaxon_pka_train_2016_03_17.pkadata");
                PKaTrainingResult pKaTrainingResult = PKaTrainingUtil.loadCorrectionData(url.openStream());
//                PKaTrainingUtil.saveCorrectionData(pKaTrainingResult, (String) null);
                plugin.setCorrectionData(pKaTrainingResult);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        plugin.setBasicpKaLowerLimit(-5);
        plugin.setAcidicpKaUpperLimit(25);
        for(PropertyMolecule mol:mols){
            if(progressReporter!=null){
                int progress = 100*mols.indexOf(mol)/mols.size();
                progressReporter.reportProgress("Progress",progress);
            }
            OEGraphMol oemol = mol.getMol();
            String molString = getMolString(oemol);
            if(molString!=null&&!molString.isEmpty()){
                try {
                    Molecule chemaxon_mol = MolImporter.importMol(molString, "sdf");
                    plugin.setMolecule(chemaxon_mol);
                    try {
                        plugin.run();
                    } catch (Exception e) {
                        e.printStackTrace();
                        continue;
//                        try {
//                            PrintWriter writer = new PrintWriter("/Users/jfeng1/pKa_error.mol");
//                            writer.println(molString);
//                            writer.close();
//                            e.printStackTrace();
//                        } catch (FileNotFoundException e1) {
//                            e1.printStackTrace();
//                        }
                    }
                    double[] acidpKas = new double[3];
                    int[] acidicIndexes = new int[3];
                    plugin.getMacropKaValues(pKaPlugin.ACIDIC,acidpKas,acidicIndexes);
                    Vector<pKa> acidic_pKas = new Vector<pKa>();
                    for(int i=0;i<3;i++){
                        int atom_id = acidicIndexes[i];
                        if(atom_id >=0){
                            pKa pKa = new pKa(atom_id,mol,new Double(acidpKas[i]).floatValue(),pKa_type.ACIDIC);
                            mol.addPKa(atom_id,pKa,PropertyMolecule.CHEMAXON_PKA);
                            acidic_pKas.add(pKa);
                        }
                    }
                    if(!acidic_pKas.isEmpty()){
                        mol.addProperty("ChemAxon Acidic pKa",""+acidic_pKas.get(0).getValue());
                        for(int i=0;i<acidic_pKas.size();i++){
                            int id = i+1;
                            String tag1 = "ChemAxon Acidic pKa " + id;
                            oechem.OESetSDData(oemol, tag1,""+acidic_pKas.get(i).getValue());
                            String tag2 = "ChemAxon Acidic pKa AtomId " + id;
                            int atomid = acidic_pKas.get(i).getAtom_id()+1;
                            oechem.OESetSDData(oemol, tag2,""+atomid);
                        }
                    }

                    Vector<pKa> basic_pKas = new Vector<pKa>();
                    double[] basicpKas = new double[3];
                    int[] basicIndexes = new int[3];
                    plugin.getMacropKaValues(pKaPlugin.BASIC,basicpKas,basicIndexes);
                    for(int i=0;i<3;i++){
                        int atom_id = basicIndexes[i];
                        if(atom_id >=0){
                            pKa pKa = new pKa(atom_id,mol,new Double(basicpKas[i]).floatValue(),pKa_type.BASIC);
                            mol.addPKa(atom_id,pKa,PropertyMolecule.CHEMAXON_PKA);
                            basic_pKas.add(pKa);
                        }
                    }
                    if(!basic_pKas.isEmpty()){
                        mol.addProperty("ChemAxon Basic pKa",""+basic_pKas.get(0).getValue());
                        for(int i=0;i<basic_pKas.size();i++){
                            int id = i+1;
                            String tag1 = "ChemAxon Basic pKa " + id;
                            oechem.OESetSDData(oemol, tag1,""+basic_pKas.get(i).getValue());
                            String tag2 = "ChemAxon Basic pKa AtomId " + id;
                            int atomid = basic_pKas.get(i).getAtom_id()+1;
                            oechem.OESetSDData(oemol, tag2,""+atomid);
                        }
                    }

                } catch (MolFormatException e) {
                    e.printStackTrace();
                } catch (PluginException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public static void calculateChemAxonLogP(Vector<PropertyMolecule> mols, ProgressReporter progressReporter){
        logPPlugin plugin = new logPPlugin();
        plugin.setDoublePrecision(2);
        plugin.setlogPMethod(LogPMethod.CHEMAXON);
        for(PropertyMolecule mol:mols){
            if(progressReporter!=null){
                int progress = 100*mols.indexOf(mol)/mols.size();
                progressReporter.reportProgress("Progress",progress);
            }
            OEGraphMol oemol = mol.getMol();
            String molString = getMolString(oemol);
            if(molString!=null&&!molString.isEmpty()){
                try {
                    Molecule chemaxon_mol = MolImporter.importMol(molString, "sdf");
                    plugin.setMolecule(chemaxon_mol);
                    plugin.run();
                    double v = plugin.getlogPTrue();
                    mol.addProperty("ChemAxon LogP", ""+v);
                } catch (MolFormatException e) {
                    e.printStackTrace();
                } catch (PluginException e) {
                    e.printStackTrace();
                } catch (ArrayIndexOutOfBoundsException e){
                    System.out.println(molString);
                    System.out.println(mol.getName());
                    e.printStackTrace();
                }
            }
        }
    }

    public static void calculateChemAxonLogDNeutral(Vector<PropertyMolecule> mols, ProgressReporter progressReporter){
        calculateChemAxonLogD(mols,7.4, progressReporter);
    }


    public static void calculateChemAxonLogD(Vector<PropertyMolecule> mols, double pH, ProgressReporter progressReporter){
        logDPlugin plugin = new logDPlugin();
        if(InSilicoToolOptions.use_trained_chemaxon_pKa) {
            try {
                URL url = new URL("http://javelin.corp.My.com:8080/insilico/chemaxon_pka_train_2016_03_17.pkadata");
                PKaTrainingResult pKaTrainingResult = PKaTrainingUtil.loadCorrectionData(url.openStream());
                PKaTrainingUtil.saveCorrectionData(pKaTrainingResult, (String) null);
                plugin.setpKaCorrectionLibrary("chemaxon_pka_train_2016_03_17");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        plugin.setDoublePrecision(2);
        plugin.setlogPMethod(LogPMethod.CHEMAXON);
        plugin.setpH(pH);
        for(PropertyMolecule mol:mols){
            if(progressReporter!=null){
                int progress = 100*mols.indexOf(mol)/mols.size();
                progressReporter.reportProgress("Running ...",progress);
            }
            OEGraphMol oemol = mol.getMol();
            String molString = getMolString(oemol);
            if(molString!=null&&!molString.isEmpty()){
                try {
                    Molecule chemaxon_mol = MolImporter.importMol(molString, "sdf");
                    plugin.setMolecule(chemaxon_mol);
                    plugin.run();
                    double v = plugin.getlogD();
                    mol.addProperty("ChemAxon LogD", ""+v);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }


    public static void generatePlasmaProteinBinding(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException{
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String ppbResult = (String)admeClient.execute("PPB",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(ppbResult);
        int idx = 0;
        for(PropertyMolecule mol:mols) {
            String moka_id = "" + idx++;
            if (obj.containsKey(moka_id)) {
                Double result = (Double) obj.get(moka_id);
                mol.addProperty("Plasma Protein Binding (Human)",result.toString());
            }
        }
    }

    public static void generateHERG(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException {
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String hERGresult = (String)admeClient.execute("hERG",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(hERGresult);
        int idx = 0;
        for(PropertyMolecule mol:mols) {
            String moka_id = "" + idx++;
            if (obj.containsKey(moka_id)) {
                JSONArray result = (JSONArray) obj.get(moka_id);
                Double value = (Double)result.get(0);
                Double lower = (Double)result.get(1);
                Double higher = (Double) result.get(2);
                MolProperty p = new MolProperty(value.toString());
                p.setRange(lower,higher);
                mol.addPropertyObj("hERG inhibition",p);
            }
        }
    }

    public static void generateRLM(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException {
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String rlmResult = (String)admeClient.execute("RLM",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(rlmResult);
        int idx = 0;
        for(PropertyMolecule mol:mols) {
            String moka_id = "" + idx++;
            if (obj.containsKey(moka_id)) {
                JSONArray result = (JSONArray) obj.get(moka_id);
                Double value = (Double)result.get(0);
                Double lower = (Double)result.get(1);
                Double higher = (Double) result.get(2);
                MolProperty p = new MolProperty(value.toString());
                p.setRange(lower,higher);
                mol.addPropertyObj("RLM Qh%",p);
            }
        }
    }

    public static void generateEffluxSVM(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException, ConnectException{
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFStringEfflux(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String effluxResult = (String)admeClient.execute("efflux_svm",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(effluxResult);
        int idx = 0;
        for(PropertyMolecule mol:mols) {
            String moka_id = "" + idx++;
            if (obj.containsKey(moka_id)) {
                Double value = (Double) obj.get(moka_id);
                MolProperty p = new MolProperty(value.toString());
                mol.addPropertyObj("Efflux Ratio(B-A/A-B)1uM (SVM)",p);
            }
        }
    }

    public static void generateEffluxKrig(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException{
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String effluxResult = (String)admeClient.execute("efflux",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(effluxResult);
        int idx = 0;
        for(PropertyMolecule mol:mols) {
            String moka_id = "" + idx++;
            if (obj.containsKey(moka_id)) {
                JSONArray result = (JSONArray) obj.get(moka_id);
                Double value = (Double)result.get(0);
                Double lower = (Double)result.get(1);
                Double higher = (Double) result.get(2);
                MolProperty p = new MolProperty(value.toString());
                p.setRange(lower,higher);
                mol.addPropertyObj("Efflux Ratio(B-A/A-B)1uM (Kriging)",p);
            }
        }
    }

    public static void toggleMDLChiralFlags(Vector<PropertyMolecule> mols, boolean on){
        for(PropertyMolecule mol:mols){
            if(mol.isSelected()) {
                oechem.OEMDLSetParity(mol.getMol(), on);
                mol.setNeedUpdate(true);
            }
        }
    }

    public static void addOrChiralFlag(Vector<PropertyMolecule> mols){
        for(PropertyMolecule mol:mols){
            if(mol.isSelected()) {
                //for group in mol.GetGroups(OEHasGroupType(OEGroupType_MDLOrStereo)):
                OEGraphMol oemol = mol.getMol();
                oechem.OEPerceiveChiral(oemol);
                oechem.OEMDLPerceiveBondStereo(oemol);
                for(OEGroupBase group:oemol.GetGroups(new OEHasGroupType(OEGroupType.MDLAndStereo))){
//                    OEAtomBaseVector v = new OEAtomBaseVector();
//                    for(OEAtomBase atom:group.GetAtoms()){
//                        v.add(atom);
//                    }
                    oemol.DeleteGroup(group);
//                    oemol.NewGroup(OEGroupType.MDLOrStereo,v);
                }
                for(OEGroupBase group:oemol.GetGroups(new OEHasGroupType(OEGroupType.MDLOrStereo))){
                    oemol.DeleteGroup(group);
                }
                for (OEAtomBase atom : oemol.GetAtoms()) {
                    if (atom.IsChiral() && atom.IsCarbon() && atom.HasStereoSpecified(OEAtomStereo.Tetrahedral)) {
                        OEAtomBaseVector v = new OEAtomBaseVector();
                        for (OEAtomBase nbr : atom.GetAtoms()) {
                            v.add(nbr);
                        }
                        int stereo = atom.GetStereo(v, OEAtomStereo.Tetrahedral);
                        System.out.println(stereo==OEAtomStereo.LeftHanded?"Left":"Right");
                        if (stereo == OEAtomStereo.LeftHanded||stereo==OEAtomStereo.RightHanded) {
                            OEAtomBaseVector orGroup = new OEAtomBaseVector();
                            orGroup.add(atom);
                            oemol.NewGroup(OEGroupType.MDLOrStereo, orGroup);
                        }
                    }
                }
                mol.setNeedUpdate(true);
            }
        }
    }

    public static void removeOrChiralFlags(Vector<PropertyMolecule> mols){
        for(PropertyMolecule mol:mols){
            if(mol.isSelected()) {
                //for group in mol.GetGroups(OEHasGroupType(OEGroupType_MDLOrStereo)):
                OEGraphMol oemol = mol.getMol();
                oechem.OEPerceiveChiral(oemol);
                oechem.OEMDLPerceiveBondStereo(oemol);
                for(OEGroupBase group:oemol.GetGroups(new OEHasGroupType(OEGroupType.MDLOrStereo))){
                    oemol.DeleteGroup(group);
                }
                mol.setNeedUpdate(true);
            }
        }
    }


    public static void predictAMES(Vector<PropertyMolecule> mols) throws MalformedURLException, XmlRpcException, ParseException {
        admeClient = getADMEClient();
        String inputSdf = convertMolVectorToSDFString(mols);
        String input = Base64.encodeBase64String(inputSdf.getBytes());
        Object[] args = new Object[]{input};
        String amesResult = (String)admeClient.execute("predictAMES",args);
        JSONParser parser = new JSONParser();
        JSONObject obj = (JSONObject)parser.parse(amesResult);
        int idx = 0;
        for(PropertyMolecule mol:mols) {
            String moka_id = "" + idx++;
            if (obj.containsKey(moka_id)) {
                String result = (String) obj.get(moka_id);
                mol.addProperty("Ames",result);
            }
        }
    }

    public static void calculateRdkitProperties(List<PropertyMolecule> propertyMolecules){
//        "hydrogen-bond acceptors",
//                "hydrogen-bond donors",
//                "molecular weight",
//                "sum of formal charges",
//                "number of rings",
//                "rotatable bonds",
//                "CLogP",
//                "2d PSA",
        NumberFormat nf = new DecimalFormat("#.##");
        for(PropertyMolecule propertyMolecule:propertyMolecules){
            RWMol rwMol = null;
            try {
                rwMol = RWMol.MolFromMolBlock(propertyMolecule.getSdfStr());
            } catch (MolSanitizeException e) {
                e.printStackTrace();
                continue;
            }
            double numHba = RDKFuncs.calcNumHBA(rwMol);
            double numHbd = RDKFuncs.calcNumHBD(rwMol);
            double mw = RDKFuncs.calcAMW(rwMol);
            double sum_of_formalchargs = RDKFuncs.getFormalCharge(rwMol);
            double num_of_rings = RDKFuncs.calcNumRings(rwMol);
            double nRots = RDKFuncs.calcNumRotatableBonds(rwMol);
            double logp = RDKFuncs.calcMolLogP(rwMol);
            double psa = RDKFuncs.calcTPSA(rwMol);
            propertyMolecule.addProperty("hydrogen-bond acceptors",nf.format(numHba));
            propertyMolecule.addProperty("hydrogen-bond donors",nf.format(numHbd));
            propertyMolecule.addProperty("molecular weight",nf.format(mw));
            propertyMolecule.addProperty("sum of formal charges",nf.format(sum_of_formalchargs));
            propertyMolecule.addProperty("number of rings",nf.format(num_of_rings));
            propertyMolecule.addProperty("rotatable bonds",nf.format(nRots));
            propertyMolecule.addProperty("CLogP",nf.format(logp));
            propertyMolecule.addProperty("2d PSA",nf.format(psa));
            rwMol.delete();
        }
    }

    public static void pickDiverseNMolecules(PropertyMolecule[] propertyMols, int n){
        if(propertyMols==null||propertyMols.length<=n||n<=0){
            return;
        }

        for(PropertyMolecule mol:propertyMols){
            mol.setIsSelected(false);
        }
        HashMap<Integer,PropertyMolecule> molHash = new HashMap<>();
        EBV_Vect fps = new EBV_Vect();
        for(PropertyMolecule mol:propertyMols){
            ROMol roMol = null;
            try {
                roMol = RWMol.MolFromMolBlock(mol.getSdfStr());
            } catch (MolSanitizeException e) {
                e.printStackTrace();
                continue;
            }
            ExplicitBitVect bitVect = RDKFuncs.getMorganFingerprintAsBitVect(roMol, 3, 2048);
            fps.add(bitVect);
            int idx = (int)(fps.size()-1);
            molHash.put(idx,mol);
        }
        if(fps.size()>n) {
            Int_Vect int_vect = RDKFuncs.pickUsingFingerprints(fps, n);
            for (int i = 0; i < int_vect.size(); i++) {
                int idx = int_vect.get(i);
                molHash.get(idx).setIsSelected(true);
            }
        }
    }

    public static void calculateOEProperty(List<PropertyMolecule> propertyMols){
        if(!oechem.OEChemIsLicensed("filter")){
            calculateRdkitProperties(propertyMols);
            return;
        }
        if(propertyMols!=null&&propertyMols.size()>0){
            OEFilter filter = new OEFilter(OEFilterType.Lead);
            oechem.OEThrow.SetLevel(OEErrorLevel.Warning);

            oeosstream ostr = new oeosstream();
            filter.SetTable(ostr);
            List<String> fields = Arrays.asList(ostr.str().split("\t"));
            ostr.clear();
            int [] propertyIdxes = new int[OE_PROPERTIES.length];
            for(int i=0;i<propertyIdxes.length;i++){
                propertyIdxes[i] = fields.indexOf(OE_PROPERTIES[i]);
            }
            for(PropertyMolecule propertyMol:propertyMols){
                OEGraphMol mol = new OEGraphMol(propertyMol.getMol());
                filter.call(mol);
                mol.delete();
                fields = Arrays.asList(ostr.str().split("\t"));
                ostr.clear(); // remove the header row from the stream
                try {
                    for(int i=0;i<propertyIdxes.length;i++){
                        String propertyValue = fields.get(propertyIdxes[i]);
                        propertyMol.addProperty(OE_PROPERTIES[i],propertyValue);
                    }
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }

    public static void calculateOEPropertySingleMol(PropertyMolecule propertyMol){
        if(!oechem.OEChemIsLicensed("filter")){
            Vector<PropertyMolecule> pmols = new Vector<>();
            pmols.add(propertyMol);
            calculateRdkitProperties(pmols);
            return;
        }
        if(!propertyMol.isEmpty()){
            OEFilter filter = new OEFilter(OEFilterType.Lead);
            oechem.OEThrow.SetLevel(OEErrorLevel.Warning);

            oeosstream ostr = new oeosstream();
            filter.SetTable(ostr);
            List<String> fields = Arrays.asList(ostr.str().split("\t"));
            ostr.clear();
            int [] propertyIdxes = new int[OE_PROPERTIES.length];
            for(int i=0;i<propertyIdxes.length;i++){
                propertyIdxes[i] = fields.indexOf(OE_PROPERTIES[i]);
            }

            OEGraphMol mol = new OEGraphMol(propertyMol.getMol());
            filter.call(mol);
            mol.delete();
            fields = Arrays.asList(ostr.str().split("\t"));
            ostr.clear(); // remove the header row from the stream
            for(int i=0;i<propertyIdxes.length;i++){
                String propertyValue = fields.get(propertyIdxes[i]);
                propertyMol.addProperty(OE_PROPERTIES[i],propertyValue);
            }
        }
    }

    public static String getSmiles(OEGraphMol mol){
        if(mol!=null){
            return oechem.OEMolToSmiles(mol);
        }
        return "";
    }

    public static String mol2Svg(PropertyMolecule propertyMolecule){
        try {
            Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(propertyMolecule.getMol());
//            MolHandler mh = new MolHandler(molecule);
//            mh.removeHydrogens();
            if(InSilicoToolOptions.show_pKa) {
                for (MolAtom atom : molecule.getAtomArray()) {
                    int atom_id = molecule.indexOf(atom);
                    if(propertyMolecule.getpKaMap(InSilicoToolOptions.pka_type).containsKey(atom_id)){
                        if(propertyMolecule.getpKaMap(InSilicoToolOptions.pka_type).get(atom_id).isAcidic()) {
                            atom.setExtraLabelColor(16711680);
                        }else {
                            atom.setExtraLabelColor(255);
                        }
                        DecimalFormat format = new DecimalFormat("#.##");
                        float value = propertyMolecule.getpKaMap(InSilicoToolOptions.pka_type).get(atom_id).getValue();
                        atom.setExtraLabel(format.format(value));
                    }
                }
            }

//            int atomId = atomDisplay.GetAtom().GetIdx();
//            if(propertyMolecule.getpKaMap(InSilicoToolOptions.pka_type).containsKey(atomId)){
//                OEFont font = atomDisplay.GetPropertyFont();
//                if(propertyMolecule.getpKaMap(InSilicoToolOptions.pka_type).get(atomId).isAcidic()){
//                    font.SetColor(oechem.getOERedOrange());
//                }else{
//                    font.SetColor(oechem.getOEGreenBlue());
//                }
//                OE2DPoint oe2DPoint = atomDisplay.GetPropertyOffset();
//                atomDisplay.SetPropertyOffset(new OE2DPoint(1.25*oe2DPoint.GetX(),1.25*oe2DPoint.GetY()));
//                atomDisplay.SetPropertyFont(font);
//            }
            /*
            MDocument doc = new MDocument(molecule);
            //Creating textbox for the name.
            MTextBox nameBox = new MTextBox();
            nameBox.setHorizontalAlignment(MTextBox.ALIGN_CENTER);
            nameBox.setVerticalAlignment(MTextBox.ALIGN_CENTER);
            //Corners of the textbox are determined. Below the molecule, centered horizontally to the centroid of the molecule.
            double molBottomEdge = doc.calcCenter().y - (molecule.calcHeight())-1;
            double molRightEdge = doc.calcCenter().x + molecule.calcWidth() / 2;
            double molLeftEdge = doc.calcCenter().x - molecule.calcWidth() / 2;
            MPoint topLeft = new MPoint(molLeftEdge, molBottomEdge-1);
            MPoint topRight = new MPoint(molRightEdge, molBottomEdge-1);
            MPoint bottomLeft = new MPoint(molLeftEdge, molBottomEdge);
            MPoint bottomRight = new MPoint(molRightEdge, molBottomEdge);
            nameBox.setPoints(new MPoint[] { bottomLeft, bottomRight, topRight, topLeft });
            //Inserting name into the textbox
            nameBox.setText(molecule.getName());
            //Adding the textbox to the document
            doc.addObject(nameBox);
            */
            byte[] svgBytes = MolExporter.exportToBinFormat(molecule,"svg");
            String svg = new String(svgBytes);
            return svg;
        } catch (IOException e) {
            e.printStackTrace();
            return "";
        }
    }

    public static byte[] mol2png(PropertyMolecule propertyMolecule, double width, double height){
        try {
            Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(propertyMolecule.getMol());
            byte[] pngBytes = MolExporter.exportToBinFormat(molecule,"png");
            return pngBytes;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }


    public static OEGraphMol getMolFromMolString(String molfile, int format) {
        oemolistream ifs = new oemolistream();
        ifs.SetFormat(format);
        OEGraphMol mol = new OEGraphMol();
        ifs.openstring(molfile);
        boolean result = oechem.OEReadMolecule(ifs, mol);
        ifs.close();
        if(result) {
            oechem.OEAssignAromaticFlags(mol);
            return mol;
        }else{
            return null;
        }
    }

    public static boolean check_library_compounds(PropertyMolecule mol){
        HashMap<String, Integer> libraryBadSmarts = InSilicoToolOptions.libraryBadSmarts;
        OESubSearch subSearch = new OESubSearch();
        for(String smarts:libraryBadSmarts.keySet()){
            subSearch.Init(smarts);
            int maxMatch = libraryBadSmarts.get(smarts);
            oechem.OEPrepareSearch(mol.getMol(),subSearch);
            if(maxMatch==0){
                if(subSearch.SingleMatch(mol.getMol())){
                    System.out.println(smarts);
                    return false;
                }
            }else{
                int numMatch  = 0;
                for(OEMatchBase ma:subSearch.Match(mol.getMol(),true)){
                    numMatch += 1;
                }
                if(numMatch > maxMatch){
                    System.out.println(smarts);
                    return false;
                }
            }
        }
        return true;
    }

    public static StructureAlert[] process_structure_alerts(){
        AlertRule[] libraryAlertRules = AlertRule.LIBRARY_ALERT_RULES;
        StringBuilder libraryRules = new StringBuilder();
        for(AlertRule alertRule:libraryAlertRules){
            libraryRules.append(alertRule.getSmarts()+"\t"+alertRule.getDescription()+"\t"+"Unwanted groups"+"\t"+alertRule.getMax()+"\n");
        }
        String[] smarts = (InSilicoToolOptions.metabolic_activation_smarts + InSilicoToolOptions.toxicophore_alert_smarts+libraryRules.toString()).split("\n");
        Vector<StructureAlert> alerts = new Vector<StructureAlert>();
        for(String s :smarts){
            if(s.startsWith("#")){
                continue;
            }
            String[] p = s.split("\t");
            if(p.length==3){
                alerts.add(new StructureAlert(p[0],p[1],p[2]));
            }else if(p.length==4){
                alerts.add(new StructureAlert(p[0],p[1],p[2],Integer.parseInt(p[3])));
            }
        }
        return alerts.toArray(new StructureAlert[alerts.size()]);
    }


    public static String getTareWeight(String barcode) throws MalformedURLException, XmlRpcException {
        XmlRpcClient admeClient = getADMEClient();
        return (String)admeClient.execute("getTareWeight",new Object[]{barcode});
    }

//    def getActiveSiteInteractionImage(receptor_mol, ligand_mol, width=1000, height=1000):
//            active_site = oebio.OEInteractionHintContainer(receptor_mol,ligand_mol)
//            oedocking.oedocking.OEAddDockingInteractions(active_site)
//            oegrapheme.oegrapheme.OEPrepareActiveSiteDepiction(active_site)
//            image = oedepict.OEImage(width,height)
//            center_frame = oedepict.OEImageFrame(image,width*0.8,height*0.8, oedepict.OE2DPoint(0.0,0.0))
//            legend_frame = oedepict.OEImageFrame(image,width*0.2,width*0.2,oedepict.OE2DPoint(width*0.8,0.0))
//            opts = oegrapheme.OE2DActiveSiteDisplayOptions(center_frame.GetWidth(),center_frame.GetHeight())
//            active_site_display = oegrapheme.OE2DActiveSiteDisplay(active_site,opts)
//            oegrapheme.oegrapheme.OERenderActiveSite(center_frame,active_site_display)
//            legend_opts = oegrapheme.OE2DActiveSiteLegendDisplayOptions(10,1)
//            oegrapheme.oegrapheme.OEDrawActiveSiteLegend(legend_frame,active_site_display,legend_opts)
//            return image



    public static void calculateCNSTempoScore(Vector<PropertyMolecule> molList){
        int [] tempo_type_list = new int[]{CNSTempo.LOGP_TYPE,CNSTempo.NUM_ROT_BONDS,CNSTempo.NUM_ARO_RINGS,CNSTempo.NUM_CHAINS,
                                           CNSTempo.NUM_H_ACCEPTORS,CNSTempo.CARBON_HETERO_RATIO,CNSTempo.NUM_BASIC_AMINES,CNSTempo.NUM_NON_CONJUGATED_RING_CARBON };
        HashMap<Integer,CNSTempo> tempoDict = new HashMap<Integer, CNSTempo>();
        for(Integer type:tempo_type_list){
            tempoDict.put(type, CNSTempo.getCNSTempoByType(type));
        }
        for(PropertyMolecule mol:molList){
            float logp = (float) mol.getProperty("Consensus LogP").getValue();
            float numRotBond = (float) mol.getProperty("rotatable bonds").getValue();
            float numAroRing = calculateNumAroRings(mol.getMol());
            float numChains = calculateNumChains(mol.getMol());
            float numHAcceptors = calculateNumHbondAcceptors(mol.getMol());
            float carbon_hetero_ratio = calculateHeteroAtomRatio(mol.getMol());
            float num_basic_amines = calculateNumBasicAmines(mol.getMol());
            float num_non_conjugated_ring_carbon = calulateNonNumConjugateNonAromaticCarbonInRing(mol.getMol());
            oechem.OESetSDData(mol.getMol(),"CNSmTEMPO_N_Aro_Rings",""+numAroRing);
            oechem.OESetSDData(mol.getMol(),"CNSmTEMPO_N_Chain_InSilico",""+numChains);
            oechem.OESetSDData(mol.getMol(),"CNSmTEMPO_Num_Acceptors_InSilico",""+numHAcceptors);
            oechem.OESetSDData(mol.getMol(),"CNSmTEMPO_carbon_hetero_ratio",""+carbon_hetero_ratio);
            oechem.OESetSDData(mol.getMol(),"CNSmTEMPO_Num_BASIC_CENTER",""+num_basic_amines);
            oechem.OESetSDData(mol.getMol(), "CNSmTEMPO_NumNonConjRingAtomgs",""+num_non_conjugated_ring_carbon);

            float[] values = new float[]{logp,numRotBond,numAroRing,numChains,numHAcceptors,carbon_hetero_ratio,num_basic_amines,num_non_conjugated_ring_carbon};
            float cnsTempo = 0.0f;
            for(Integer type:tempo_type_list){
                cnsTempo += tempoDict.get(type).getPenaltyCoeff(values[type]);
            }
            if(cnsTempo>10){
                cnsTempo = 10;
            }
            cnsTempo = (10-cnsTempo)/10*6.0f;
            String value = String.format("%5.2f", cnsTempo);
            mol.addPropertyObj("CNS mTEMPO", new CNSTempoObj(value));
        }

    }

    public static float calculateHeteroAtomRatio(OEGraphMol mol){
        if(mol==null){
            return 0;
        }
        float carbon_count = 0.0f;
        float non_carbon_count = 0.0f;
        for(OEAtomBase atm:mol.GetAtoms()){
            if(atm.IsHydrogen()){
                continue;
            }
            if(atm.IsCarbon()){
                carbon_count+=1;
            }else{
                non_carbon_count += 1;
            }
        }
        if(non_carbon_count==0){
            non_carbon_count = 0.01f;
        }
        return carbon_count/non_carbon_count;
    }

    public static float calulateNonNumConjugateNonAromaticCarbonInRing(OEGraphMol mol){
        if(mol==null){
            return 0.0f;
        }
        oechem.OEFindRingAtomsAndBonds(mol);
        oechem.OEAssignAromaticFlags(mol);
        String smarts = "[C;r;!a]=[C;r;!a]-[A]=[A]";
        OESubSearch subSearch = new OESubSearch();
        subSearch.Init(smarts);
        oechem.OEPrepareSearch(mol,subSearch);
        HashSet<Integer> atomids = new HashSet<Integer>();
        for(OEMatchBase match:subSearch.Match(mol)){
            for(OEMatchPairAtom ma:match.GetAtoms()){
                OEAtomBase atm = ma.getTarget();
                if(atm.IsCarbon()&&atm.IsInRing()&&!atm.IsAromatic()){
                    atomids.add(atm.GetIdx());
                }
            }
        }
        int numConjugatedCarbonAtom =  atomids.size();
        int total_carbon_atom_in_nonAroRings = 0;
        for(OEAtomBase atm:mol.GetAtoms()){
            if(!atm.IsAromatic()&&atm.IsCarbon()&&atm.IsInRing()){
                total_carbon_atom_in_nonAroRings += 1;
            }
        }
        return total_carbon_atom_in_nonAroRings-numConjugatedCarbonAtom;
    }

    public static int calculateNumBasicAmines(OEGraphMol mol){
        if(mol==null){
            return 0;
        }
        String smarts = "[#7;!$([#7]=,:,#*);!$([#7]-,:*=,:*);!$(N-N)]";
        OESubSearch subSearch = new OESubSearch();
        subSearch.Init(smarts);
        oechem.OEPrepareSearch(mol,subSearch);
        HashSet<Integer> atomids = new HashSet<Integer>();
        for(OEMatchBase match:subSearch.Match(mol)){
            for(OEMatchPairAtom ma:match.GetAtoms()){
                OEAtomBase atm = ma.getTarget();
                atomids.add(atm.GetIdx());
            }
        }
        return atomids.size();

//        Vector<PropertyMolecule> molList = new Vector<PropertyMolecule>();
//        PropertyMolecule pmol = new PropertyMolecule(mol);
//        molList.add(pmol);
//        ChemFunc.calculateChemAxon_pKa(molList,null);
//        HashMap<Integer,pKa> pKaDict = pmol.getpKaMap(PropertyMolecule.CHEMAXON_PKA);
//        int nBasicAtoms = 0;
//        for(OEAtomBase atom:pmol.getMol().GetAtoms()){
//            if(pKaDict.containsKey(atom.GetIdx())){
//                pKa pKa = pKaDict.get(atom.GetIdx());
//                if(pKa.getType()== pKa_type.BASIC&&pKa.getValue()>7.0){
//                    nBasicAtoms += 1;
//                }
//            }
//        }
//        return nBasicAtoms;
    }

    public static int calculateNumChains(OEGraphMol oemol){
        if(oemol==null){
            return 0;
        }
        int[] parts = new int[oemol.GetMaxAtomIdx()];
        oechem.OEFindRingAtomsAndBonds(oemol);
        int nRings = oechem.OEDetermineRingSystems(oemol, parts);
        Vector<Integer> chainAtomList = new Vector<Integer>();
        for(int i=0;i<parts.length;i++){
            OEAtomBase atm = oemol.GetAtom(new OEHasAtomIdx(i));
            if(atm==null||atm.IsHydrogen()){
                continue;
            }
            if(parts[i]==0) {
                chainAtomList.add(i);
            }
        }
        int nChains = 0;
        for(int atom_id:chainAtomList){
            OEAtomBase atm = oemol.GetAtom(new OEHasAtomIdx(atom_id));
            if(atm.GetHvyDegree()==1){
                nChains +=1;
            }
        }
        if(nRings==0){
            nChains -=1;
        }

        OEAtomBondSetIter oeAtomBondSets = oemedchem.OEGetBemisMurcko(oemol);
        for(OEAtomBondSet a:oeAtomBondSets){
            String roleName = null;
            for(OERole r:a.GetRoles()){
                roleName = r.GetName();
                break; //Only one role for OEGetBemisMurcko
            }
            if(roleName!=null&&roleName.equals("Linker")){
                OEGraphMol frag = new OEGraphMol();
                oechem.OESubsetMol(frag, oemol, a, true);
                int[] frag_parts = new int[frag.GetMaxAtomIdx()];
                int numFrags = oechem.OEDetermineComponents(frag,frag_parts);
                nChains += numFrags;
            }
        }
        for(OEBondBase bnd:oemol.GetBonds()){
            if(bnd.GetBgn().IsInRing()&&bnd.GetEnd().IsInRing()&&!bnd.IsInRing()){
                nChains+= 1;
            }
        }
        return nChains;
    }

    public static int calculateNumAroRings(OEGraphMol oemol){
        if(oemol==null){
            return 0;
        }
        int[] parts = new int[oemol.GetMaxAtomIdx()];
        oechem.OEFindRingAtomsAndBonds(oemol);
        oechem.OEAssignAromaticFlags(oemol);
        oechem.OEDetermineAromaticRingSystems(oemol, parts);
        HashMap<Integer,Vector<Integer>> ringMap = new HashMap<Integer, Vector<Integer>>();
        for(int i=0;i<parts.length;i++){
            int ring_idx = parts[i];
            if(ring_idx==0){
                continue;
            }
            if(!ringMap.containsKey(ring_idx)){
                ringMap.put(ring_idx,new Vector<Integer>());
            }
            ringMap.get(ring_idx).add(i);
        }
        int n_aro_rings = 0;
        for(Integer ringIdx:ringMap.keySet()){
            Vector<Integer> atomIdxList = ringMap.get(ringIdx);
            Vector<OEBondBase> bondList = new Vector<OEBondBase>();
            for(Integer atomId:atomIdxList){
                OEAtomBase atm = oemol.GetAtom(new OEHasAtomIdx(atomId));
                for(OEBondBase bond:atm.GetBonds()){
                    if(bondList.contains(bond)){
                        continue;
                    }else{
                        if(atomIdxList.contains(bond.GetBgnIdx())&&atomIdxList.contains(bond.GetEndIdx())){
                            bondList.add(bond);
                        }
                    }
                }
            }
            n_aro_rings += bondList.size()-atomIdxList.size()+1;
        }
        return n_aro_rings;
    }

    public static void testPMI(){
        oemolistream ifs = new oemolistream();
        ifs.open("/Users/jfeng1/00Demo/rod.sdf");
        OEGraphMol mol = new OEGraphMol();
        while(oechem.OEReadMolecule(ifs,mol)){
            double[] pmi = ChemFunc.calculatePMI(mol);
            System.out.println(String.format("%f %f",pmi[0],pmi[1]));
        }
        ifs.close();
    }

    public static int calculateNumHbondAcceptors(OEGraphMol mol){
        if(mol==null){
            return 0;
        }
        String hba_smarts = "[$([#7,#8]);!$([nX3,NX4,#7v5]);!$([NX3][c,C&X3,C&X2])]";
        OESubSearch subSearch = new OESubSearch();
        subSearch.Init(hba_smarts);
        oechem.OEPrepareSearch(mol,subSearch);
        HashSet<Integer> atomids = new HashSet<Integer>();
        for(OEMatchBase match:subSearch.Match(mol)){
            for(OEMatchPairAtom ma:match.GetAtoms()){
                OEAtomBase atm = ma.getTarget();
                    atomids.add(atm.GetIdx());
            }
        }
        return atomids.size();
    }

    public static Vector<OEGraphMol> splitMol(OEGraphMol mol){
        Vector<OEGraphMol> partMols = new Vector<>();
        int[] parts = new int[mol.GetMaxAtomIdx()];
        int pcount = oechem.OEDetermineComponents(mol,parts);
        OEPartPredAtom pred = new OEPartPredAtom(parts);
        for (int i = 1; i <= pcount; ++i) {
            pred.SelectPart(i);
            OEGraphMol partmol = new OEGraphMol();
            oechem.OESubsetMol(partmol, mol, pred, true);
            partMols.add(partmol);
        }
        return partMols;
    }

    public static void main(String[] args) {
        try {
            testPost();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (ParseException e) {
            e.printStackTrace();
        }
    }

    private void test_msm(){
        //todo:
        //todo:
        //todo:

    }

    private static void test_formatCY() {
        String a = "CY-2009775\n" +
                "CY-0000121\n" +
                "CY-0000076\n" +
                "CY-0009459\n" +
                "CY-0000255\n" +
                "CY-0009608\n" +
                "CY-0009686\n" +
                "CY-0010758\n" +
                "CY-0000181\n" +
                "CY-0010893\n" +
                "CY-10894\n" +
                "CY-10877\n" +
                "CY-13690\n" +
                "CY-11055\n" +
                "CY-151\n" +
                "CY-10278\n" +
                "CY-10742\n" +
                "CY-13702\n" +
                "CY-12458\n" +
                "CY-11058\n" +
                "CY-9767\n" +
                "CY-9453\n" +
                "CY-187\n" +
                "CY-299\n" +
                "CY-10879\n" +
                "CY-9791\n" +
                "CY-11056\n" +
                "CY-9759\n" +
                "CY-14201\n" +
                "CY-11191\n" +
                "CY-9798\n" +
                "CY-11041\n" +
                "CY-14202\n" +
                "CY-9782\n" +
                "CY-9799\n" +
                "CY-9785\n" +
                "CY-9789\n" +
                "CY-9764\n" +
                "CY-11044\n" +
                "CY-13702\n" +
                "CY-48\n" +
                "CY-15564\n" +
                "CY-48\n" +
                "CY-48\n" +
                "CY-0055180\n" +
                "CY-0052973\n" +
                "CY-71083\n" +
                "CY-71085\n" +
                "CY-66062\n" +
                "CY-9775\n" +
                "CY-14402\n" +
                "CY-66066\n" +
                "CY-82\n" +
                "CY-151\n" +
                "CY-51\n" +
                "CY-9681\n" +
                "CY-9719\n" +
                "CY-9818\n" +
                "CY-9738\n" +
                "CY-9945\n" +
                "CY-9742\n" +
                "CY-9744\n" +
                "CY-9740\n" +
                "CY-9777\n" +
                "CY-9778\n" +
                "CY-10762\n" +
                "CY-10758\n" +
                "CY-14248\n" +
                "CY-9778\n" +
                "CY-9738\n" +
                "CY-9738\n" +
                "CY-9738\n" +
                "CY-9738\n" +
                "CY-0011102\n" +
                "CY-0076683\n" +
                "CY-16502\n" +
                "CY-11029\n" +
                "CY-10636\n" +
                "CY-10338\n" +
                "CY-10735\n" +
                "CY-10176\n" +
                "CY-10145\n" +
                "CY-10377\n" +
                "CY-10356\n" +
                "CY-10618\n" +
                "CY-10328\n" +
                "CY-10482\n" +
                "CY-10424\n" +
                "CY-10479\n" +
                "CY-10249\n" +
                "CY-10059\n" +
                "4451\n";
        String[] list = a.split("\n");
        for(String tmp:list) {
            System.out.println(formatCYNumber(tmp));
        }
    }

    private static void testPost() throws IOException, ParseException {
        URL obj = new URL("http://herg.cellarity.int/predict");
        HttpURLConnection con = (HttpURLConnection) obj.openConnection();
        con.setRequestProperty("accept","application/json");
        con.setRequestProperty("Content-Type", "application/json");
        con.setRequestMethod("POST");
        JSONParser parser = new JSONParser();
        JSONArray smilesList = new JSONArray();

        JSONObject smilesObj1 = new JSONObject();
        smilesObj1.put("id",0);
        smilesObj1.put("smiles","CCCCCCCCCN");
        smilesList.add(smilesObj1);

        JSONObject smilesObj2 = new JSONObject();
        smilesObj2.put("id",1);
        smilesObj2.put("smiles","CCCCCC(=O)O");
        smilesList.add(smilesObj2);

        String parmas = smilesList.toJSONString();

        // For POST only - START
        con.setDoOutput(true);
        OutputStream os = con.getOutputStream();
        os.write(parmas.getBytes());
        os.flush();
        os.close();
        // For POST only - END

        int responseCode = con.getResponseCode();
        System.out.println("POST Response Code :: " + responseCode);

        if (responseCode == HttpURLConnection.HTTP_OK) { //success
            BufferedReader in = new BufferedReader(new InputStreamReader(con.getInputStream()));
            String inputLine;
            StringBuffer response = new StringBuffer();

            while ((inputLine = in.readLine()) != null) {
                response.append(inputLine);
            }
            in.close();

            // print result
            JSONObject result = (JSONObject) parser.parse(response.toString());
            System.out.println(result.toString());
            for(Object key : result.keySet()){
                int compound_id = Integer.parseInt((String) key);
                JSONObject molResult = (JSONObject) result.get(key);
                for(Object propertyName: molResult.keySet()){
                    System.out.println(propertyName);
                    System.out.println(molResult.get(propertyName));
                }
            }
        } else {
            System.out.println("POST request did not work.");
        }
    }


    protected static void test_canonicalization() {
        String mol_string = "VIR-0089945\n" +
                "  Mrv2003 08242005212D          \n" +
                "\n" +
                " 33 38  0  0  0  0            999 V2000\n" +
                "   -5.2826   -0.3415    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -4.4670   -0.4297    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -3.9812    0.2383    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -4.3139    0.9944    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -5.1350    1.0749    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -5.6193    0.4070    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -3.1616    0.1559    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -2.8268   -0.6011    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -3.3125   -1.2693    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -4.1333   -1.1804    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -6.4353    0.4909    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -6.7768    1.2413    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -7.5971    1.3255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -8.0758    0.6593    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -7.7425   -0.0903    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -6.9223   -0.1745    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -4.6166   -1.8431    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -2.0111   -0.6880    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.4606   -0.0799    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.6953    0.4120    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.2097    1.0803    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.2097   -0.2563    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.4240   -0.0010    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    1.4240    0.8251    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.7120    1.2418    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.8251    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.7098   -0.4112    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    0.7120   -0.4083    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    3.5156    0.4121    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "    2.4632    1.8604    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -0.7963   -1.2307    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "   -1.6004   -1.3995    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n" +
                "  1  2  2  0  0  0  0\n" +
                "  1  6  1  0  0  0  0\n" +
                "  2  3  1  0  0  0  0\n" +
                "  2 10  1  0  0  0  0\n" +
                "  3  4  2  0  0  0  0\n" +
                "  3  7  1  0  0  0  0\n" +
                "  4  5  1  0  0  0  0\n" +
                "  5  6  2  0  0  0  0\n" +
                "  6 11  1  0  0  0  0\n" +
                "  7  8  1  0  0  0  0\n" +
                "  8  9  2  0  0  0  0\n" +
                "  8 18  1  0  0  0  0\n" +
                "  9 10  1  0  0  0  0\n" +
                " 10 17  2  0  0  0  0\n" +
                " 11 12  1  0  0  0  0\n" +
                " 11 16  1  0  0  0  0\n" +
                " 12 13  1  0  0  0  0\n" +
                " 13 14  1  0  0  0  0\n" +
                " 14 15  1  0  0  0  0\n" +
                " 15 16  1  0  0  0  0\n" +
                " 18 19  2  0  0  0  0\n" +
                " 18 33  1  0  0  0  0\n" +
                " 19 28  1  0  0  0  0\n" +
                " 20 21  1  0  0  0  0\n" +
                " 20 22  1  0  0  0  0\n" +
                " 20 30  2  0  0  0  0\n" +
                " 21 24  1  0  0  0  0\n" +
                " 21 31  1  0  0  0  0\n" +
                " 22 23  1  0  0  0  0\n" +
                " 23 24  2  0  0  0  0\n" +
                " 23 29  1  0  0  0  0\n" +
                " 24 25  1  0  0  0  0\n" +
                " 25 26  2  0  0  0  0\n" +
                " 26 27  1  0  0  0  0\n" +
                " 27 28  1  0  0  0  0\n" +
                " 27 29  2  0  0  0  0\n" +
                " 28 32  1  0  0  0  0\n" +
                " 32 33  2  0  0  0  0\n" +
                "M  END\n" +
                "\n";
        try {
            String new_mol_string = ChemFunc.getCanonicalizedStructure(mol_string);
            System.out.println(new_mol_string);
        } catch (MalformedURLException e) {
            e.printStackTrace();
        } catch (XmlRpcException e) {
            e.printStackTrace();
        }
    }

    private static void test_cnsmpo() {
        oemolostream ofs = new oemolostream();
        ofs.open("/Users/jfeng1/Datasets/CNS_TEMPO/cns_tempo_tmp.sdf");
        oemolistream ifs = new oemolistream();
        ifs.open("/Users/jfeng1/Datasets/CNS_TEMPO/cns_tempo_all.sdf");
        OEGraphMol mol = new OEGraphMol();
        while(oechem.OEReadMolecule(ifs,mol)){
            float numAroRing = calculateNumAroRings(mol);
            float numChains = calculateNumChains(mol);
            float numHAcceptors = calculateNumHbondAcceptors(mol);
            float carbon_hetero_ratio = calculateHeteroAtomRatio(mol);
            float num_basic_amines = calculateNumBasicAmines(mol);
            float num_non_conjugated_ring_carbon = calulateNonNumConjugateNonAromaticCarbonInRing(mol);
            oechem.OESetSDData(mol,"N_Aro_Rings_InSilico",""+numAroRing);
            oechem.OESetSDData(mol,"N_Chain_InSilico",""+numChains);
            oechem.OESetSDData(mol,"Num_Acceptors_InSilico",""+numHAcceptors);
            oechem.OESetSDData(mol,"carbon_hetero_ratio",""+carbon_hetero_ratio);
            oechem.OESetSDData(mol,"Num_BASIC_CENTER_ISWT",""+num_basic_amines);
            oechem.OESetSDData(mol, "NumNonConjRingAtomgs",""+num_non_conjugated_ring_carbon);
            oechem.OEWriteMolecule(ofs,mol);
        }
        ofs.close();
        ifs.close();
    }
}
