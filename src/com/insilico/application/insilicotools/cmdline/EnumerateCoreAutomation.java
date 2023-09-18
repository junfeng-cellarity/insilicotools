package com.insilico.application.insilicotools.cmdline;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.util.Enumerator;
import openeye.oechem.*;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.util.HashMap;
import java.util.Vector;


/**
 * Created by jfeng1 on 2/24/16.
 */
public class EnumerateCoreAutomation {
    public static String coreDirectory = "/home/jfeng/datasets/VirtualIdea/Pharmaron/";
    public static String reagentDirectory = "/home/jfeng/datasets/BuildingBlocks/Pharmaron/";
    HashMap<String,String> reagentDict = new HashMap<String, String>();
    public static HashMap<String,String> reagentListLibrary = new HashMap<String, String>(){
        {
            put("BocDiamines", reagentDirectory+"Boc-Diamine.sdf");
            put("BoronicAcidsAndEsters", reagentDirectory+"BoronicAcids.sdf");
            put("ArylAmines", reagentDirectory+"ArylAmine.sdf");
            put("Halides", reagentDirectory+"ArylBrI.sdf");
        }
    };

    public static HashMap<String,String> smirksLibrary = new HashMap<String, String>(){
        {
            put("Sky-P-ML-22_BocDiamines_ArylAmines","[#6:2]-[#7;X3;H2,H1;!$(NC=O):1].[#6;a:3]-[#7;H2X3;!$(NC=O):4]>>[#6:3]-[#7:4]-[#6](=O)-[#6]-1=[#6](F)-[#6]=[#6](-[#6](-[#8])=[#6]-1)-[#6]-1=[#7]-[#7]=[#6](-[#7:1]-[#6:2])-[#6]=[#6]-1");
        }
    };

    public static Vector<PropertyMolecule> getReagents(String file, Vector<String> excludedList){
        Vector<PropertyMolecule> molecules = new Vector<PropertyMolecule>();
        try {
            oemolistream ifs  = new oemolistream();
            ifs.open(file);
            OEGraphMol mol = new OEGraphMol();
            while(oechem.OEReadMolecule(ifs,mol)){
                PropertyMolecule pmol = new PropertyMolecule(new OEGraphMol(mol));
                if(excludedList==null){
                    molecules.add(pmol);
                }else{
                    if(!excludedList.contains(pmol.getName())){
                        molecules.add(pmol);
                    }
                }
            }
            ifs.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return molecules;
    }

    public static double getMaxPSA(Vector<PropertyMolecule> propertyMolecules){
        double maxPSA = 0.0;
        for(PropertyMolecule m:propertyMolecules){
            if(m.getPSA()>maxPSA){
                maxPSA = m.getPSA();
            }
        }
        return maxPSA;
    }

    public static double getMaxMW(Vector<PropertyMolecule> propertyMolecules){
        double maxMW = 0.0;
        for(PropertyMolecule m:propertyMolecules){
            if(m.getMW()>maxMW){
                maxMW = m.getMW();
            }
        }
        return maxMW;
    }

    public static void main(String[] args1) {
        oechem.OEThrow.SetLevel(OEErrorLevel.Error);
        LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
        Logger rootLogger = loggerContext.getLogger("chemaxon");
        rootLogger.setLevel(Level.OFF);

        int n = 0;
        for(String key:smirksLibrary.keySet()){
            String[] args = key.split("_");
            String core = args[0];
            String reagentType1 = args[1];
            String reagentType2 = args[2];
            oemolostream ofs = new oemolostream();
            String fname = String.format("%s%s/%s_enumerate.sdf",coreDirectory,core,core);
            ofs.open(fname);
            boolean deprotect = false;
            if(reagentType1.equals("BocDiamines")||reagentType2.equals("BocDiamines")){
                deprotect = true;
            }
            String smirks = smirksLibrary.get(key);
            Enumerator enumerator = new Enumerator();
            enumerator.setSmirks(smirks);
            Vector<PropertyMolecule> reagentList1 = getReagents(reagentListLibrary.get(reagentType1), new Vector<>());
            Vector<PropertyMolecule> reagentList2 = getReagents(reagentListLibrary.get(reagentType2), new Vector<>());
            enumerator.addReagent(reagentList1,0);
            enumerator.addReagent(reagentList2,1);
            try {
                enumerator.enumerateToFile(ofs, deprotect, new ProgressReporter() {
                    @Override
                    public void reportProgress(String note, int progress) {
                        System.out.println(String.format("%s %d compounds",note,progress));
                    }
                });
                ofs.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}
