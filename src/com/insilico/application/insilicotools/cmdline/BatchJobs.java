package com.insilico.application.insilicotools.cmdline;

import chemaxon.license.LicenseManager;
import chemaxon.license.LicenseProcessingException;
import com.insilico.application.insilicotools.InSilicoToolOptions;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.data.SerializableMol;
import com.insilico.application.insilicotools.data.StructureAlert;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.util.OEChemWebLicenseInstaller;
import com.google.common.io.Files;
import openeye.oechem.*;
import org.ehcache.Cache;
import org.ehcache.CacheManager;
import org.ehcache.config.builders.CacheConfigurationBuilder;
import org.ehcache.config.builders.CacheManagerBuilder;
import org.ehcache.config.builders.ResourcePoolsBuilder;
import org.ehcache.config.units.MemoryUnit;

import javax.swing.*;
import java.io.IOException;
import java.util.Arrays;
import java.util.Vector;
import java.util.concurrent.*;

/**
 * Created by jfeng1 on 10/12/17.
 */
public class BatchJobs {
    public static int NUM_CORES = Runtime.getRuntime().availableProcessors();
//    static String[] PROPERTIES_TO_USE = {"CLogP","molecular weight","2d PSA","CNS mTEMPO",
//            "CNS MPO","ChemAxon Acidic pKa","ChemAxon Basic pKa","Structure Alerts",};
    static String[] PROPERTIES_TO_USE = {"CNS MPO","CNS mTEMPO","ChemAxon LogD","ChemAxon Acidic pKa","ChemAxon Basic pKa","MoKa Basic pKa","MoKa Acidic pKa","MoKa LogD","cHLM","cHPPB","cMDR1","cRLM","cRPPB","cSolubility"};
    private static ExecutorService executor = Executors.newFixedThreadPool(NUM_CORES);
    final static int batch_size = 100;
    public static void calculateProperties_all(String filename, String outputfile, final Vector<String> selectedProperties, ProgressReporter progressReporter) throws IOException {
        if(progressReporter!=null){
            progressReporter.reportProgress("Running",-1);
        }
        Vector<Future<?>> futures = new Vector<Future<?>>();
        oemolistream ifs = new oemolistream(filename);
        OEMolDatabase db = new OEMolDatabase(ifs);
        int max_id = db.GetMaxMolIdx();
        int id = 0;
        CacheManager cacheManager = CacheManagerBuilder.newCacheManagerBuilder().with(CacheManagerBuilder.persistence(Files.createTempDir()))
                .withCache("preConfigured",
                        CacheConfigurationBuilder.newCacheConfigurationBuilder(Integer.class, SerializableMol.class,
                                ResourcePoolsBuilder.heap(2*max_id).offheap(1,MemoryUnit.GB).disk(10,MemoryUnit.GB))
                                .build())
                .build(true);

        final Cache<Integer, SerializableMol> myCache= cacheManager.getCache("preConfigured", Integer.class, SerializableMol.class);


//        final Cache<Integer, String> myCache = cacheManager.createCache("myCache",
//                CacheConfigurationBuilder.newCacheConfigurationBuilder(Integer.class, String.class,
//                        ResourcePoolsBuilder.heap(100)).build());
        final Vector<Integer> molIds = new Vector<>();
        OEGraphMol mol = new OEGraphMol();
        while(oechem.OEReadMolecule(ifs,mol)){
            myCache.put(id,new PropertyMolecule(mol).getSerializableMol2D());
            molIds.add(id);

            if(id!=0&&(id%batch_size==0||id==max_id-1)){
                final Vector<Integer> pmols = new Vector<>();
                final int future_id = futures.size();
                futures.add(executor.submit(new Runnable() {
                    @Override
                    public void run() {
                        try {
                            calculateProperties(selectedProperties, myCache, future_id, batch_size);
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }));
            }
            id++;
        }
        ifs.close();

        try {
            while(true){
                int flag = 0;
                for(int i=0;i<futures.size();i++){
                    Future<?> f = futures.get(i);
                    if(!f.isDone()||myCache.get(i)==null){
                        flag = 1;
                        try {
                            Thread.sleep(10);
                        } catch (InterruptedException e) {
                            e.printStackTrace();
                        }
                    }
                }
                if(flag==0) {
                    break;
                }
            }
            oemolostream ofs = new oemolostream();
            ofs.open(outputfile);
            for(int molid:molIds){
                SerializableMol p = myCache.get(molid);
                oechem.OEWriteMolecule(ofs,p.getOEMol());
            }
            ofs.close();
            executor.shutdown();
            if(progressReporter!=null){
                progressReporter.reportProgress("finished",100);
            }
            System.out.println("Future finished ...");
            cacheManager.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static void calculateProperties(Vector<String> selectedProperties, Cache<Integer, SerializableMol> myCache, int future_id, int trunk_size) throws Exception{
        System.out.println("calculating ...");
        final long startTime = System.currentTimeMillis();
        Vector<PropertyMolecule> molecules = new Vector<>();
        for(int i=0;i<trunk_size;i++){
            int id = future_id * trunk_size + i;
            SerializableMol serializableMol = myCache.get(id);
            if(serializableMol !=null) {
                molecules.add(new PropertyMolecule(serializableMol.getOEMol()));
            }
        }
        if(molecules.isEmpty()){
            return;
        }

        if(selectedProperties.contains("cHLM")||selectedProperties.contains("cHPPB")||selectedProperties.contains("cMDR1")||selectedProperties.contains("cRLM")||selectedProperties.contains("cRPPB")||selectedProperties.contains("cSolubility")){
            ChemFunc.generateInSilicoADMEproperties(molecules);
        }

        ChemFunc.calculateOEProperty(molecules);
        if(selectedProperties.contains("Consensus LogP")||selectedProperties.contains("CNS mTEMPO")||selectedProperties.contains("ChemAxon LogP")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")){
            ChemFunc.calculateChemAxonLogP(molecules,null);
        }
        if(selectedProperties.contains("Consensus LogP")||selectedProperties.contains("CNS mTEMPO")){
            for(PropertyMolecule mol:molecules){
                mol.addConsensusLogP();
            }
        }
        if(selectedProperties.contains("CNS mTEMPO")){
            ChemFunc.calculateCNSTempoScore(molecules);
        }
        if(selectedProperties.contains("ChemAxon LogD")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")){
            ChemFunc.calculateChemAxonLogDNeutral(molecules,null);
        }
        if(selectedProperties.contains("ChemAxon Acidic pKa")||selectedProperties.contains("ChemAxon Basic pKa")||selectedProperties.contains("CNS MPO")||selectedProperties.contains("CNS mTEMPO")){
            ChemFunc.calculateChemAxon_pKa(molecules,null);
            InSilicoToolOptions.pka_type = PropertyMolecule.CHEMAXON_PKA;
        }
        if(selectedProperties.contains("CNS MPO")){
            ChemFunc.generateCNSMPODescriptors(molecules, false);
        }

        if(selectedProperties.contains("No. Aromatic Rings")){
            ChemFunc.calculateNoAroRings(molecules);
        }
        if(selectedProperties.contains("Plasma Protein Binding (Human)")){
            ChemFunc.generatePlasmaProteinBinding(molecules);
        }
        if(selectedProperties.contains("Structure Alerts")){
            StructureAlert[] alerts = ChemFunc.process_structure_alerts();
            for(PropertyMolecule mol:molecules){
                mol.addStructureAlerts(alerts);
            }
        }
        if(selectedProperties.contains("hERG inhibition")){
            ChemFunc.generateHERG(molecules);
        }

        if(selectedProperties.contains("RLM Qh%")){
            ChemFunc.generateRLM(molecules);
        }

        if(selectedProperties.contains("Efflux Ratio(B-A/A-B)1uM")){
            ChemFunc.generateEffluxKrig(molecules);
            ChemFunc.generateEffluxSVM(molecules);
        }


        if(selectedProperties.contains("Human VDss(L/kg)")){
            ChemFunc.generateVDss(molecules);
        }

        if(selectedProperties.contains("Ames")){
            ChemFunc.predictAMES(molecules);
        }

        if(selectedProperties.contains("Clearance Route")){
            ChemFunc.predictClearanceMechanism(molecules);
            selectedProperties.add("PC_1");
            selectedProperties.add("PC_2");
        }

        //"Metabolic Clearance","Renal Clearance"
        if(selectedProperties.contains("Human Renal Clearance(mL/min/kg)")){
            ChemFunc.predictClearance(molecules,false);
        }

        if(selectedProperties.contains("Human Metabolic Clearance(mL/min/kg)")){
            ChemFunc.predictClearance(molecules,true);
        }

        if(selectedProperties.contains("Principal Moment of Inertia (NPR1,NPR2)")){
            ChemFunc.calculatePMIBatch(molecules);
            selectedProperties.remove("Principal Moment of Inertia (NPR1,NPR2)");
            selectedProperties.add("NPR1");
            selectedProperties.add("NPR2");
        }

        //"Volume3D", "SurfaceArea3D","PolarSurfaceArea3D", "DipoleMoment3D"
        if(selectedProperties.contains("Volume3D")||selectedProperties.contains("SurfaceArea3D")
                ||selectedProperties.contains("PolarSurfaceArea3D")||selectedProperties.contains("DipoleMoment3D")){
            for(PropertyMolecule m:molecules) {
                OEChemFunc.getInstance().calculate3DSurfaceVolumeDipole(m);
            }
        }
        final long duration = System.currentTimeMillis() - startTime;
        for(int i=0;i<trunk_size;i++){
            int id = future_id * trunk_size + i;
            if(myCache.containsKey(id)&&(i<molecules.size())) {
                try {
                    myCache.put(id, molecules.get(i).getSerializableMol2D());
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        molecules.clear();
    }


    public static void main(String[] args) {
        if(args==null||args.length!=2){
            System.out.println("Usage:this_program input.sdf output.sdf");
            System.exit(1);
        }
        System.out.println(args[0]+" "+args[1]);
        try {
            OEChemWebLicenseInstaller.loadOELicenseFromWeb();
//                    LicenseManager.setLicenseFile("/Users/jfeng1/.chemaxon/license.cxl");
            LicenseManager.setLicenseFile("http://10.74.2.128:8080/insilico_tools/license.cxl");
        } catch (IOException e) {
            JOptionPane.showMessageDialog(null,e.getMessage());
            return;
        } catch (LicenseProcessingException e) {
            JOptionPane.showMessageDialog(null,e.getMessage());
            return;
        }

        final long startTime = System.currentTimeMillis();
        final Vector<String> selectedProperties = new Vector<>();
        selectedProperties.addAll(Arrays.asList(PROPERTIES_TO_USE));

        try {
            calculateProperties_all(args[0], args[1], selectedProperties, new ProgressReporter() {
                @Override
                public void reportProgress(String note, int progress) {
                    System.out.println(note+" "+progress);
                }
            });
        } catch (IOException e) {
            e.printStackTrace();
        }
        final long duration = System.currentTimeMillis() - startTime;
        System.out.println("Running time:"+duration/1000+ " secs.");
    }
}
