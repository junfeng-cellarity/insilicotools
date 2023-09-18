package com.insilico.application.insilicotools.util;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oegraphsim.*;

import java.util.Random;
import java.util.Vector;

/**
 * Created by jfeng1 on 2/20/16.
 */
public class StatFunc {
    public static Vector<PropertyMolecule> filterByMwAndPSA(Vector<PropertyMolecule> molecules, double max_mw, double max_psa){
        Vector<PropertyMolecule> filteredMolecules = new Vector<PropertyMolecule>();
        for(PropertyMolecule m:molecules){
            if(m.getMW()<max_mw&&m.getPSA()<max_psa){
                filteredMolecules.add(m);
            }
        }
        return filteredMolecules;
    }

    public static void shuffle_array(int[] ar){
        Random rnd = new Random(System.currentTimeMillis());
        for (int i = ar.length - 1; i > 0; i--)
        {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            int a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }

    public static PropertyMolecule getMolFromDifferentCluster(Vector<PropertyMolecule> currentMols, Vector<PropertyMolecule> mols, double psa, double mw){
        String clusterTag = "KMeanCluster_MACCS166";
        Vector<String> clustersList = new Vector<String>();
        for(PropertyMolecule p:currentMols){
            clustersList.add(oechem.OEGetSDData(p.getMol(),clusterTag));
        }
        int[] ar = new int[mols.size()];
        for(int i=0;i<mols.size();i++){
            ar[i] = i;
        }
        shuffle_array(ar);

        for(int i=0;i<ar.length;i++){
            PropertyMolecule p = mols.get(ar[i]);
            if(psa>0){
                double r_psa = Double.parseDouble(oechem.OEGetSDData(p.getMol(),"R_PSA"));
                if(r_psa>psa){
                    continue;
                }
            }
            if(mw>0){
                double r_mw = Double.parseDouble(oechem.OEGetSDData(p.getMol(),"R_MW"));
                if(r_mw>mw){
                    continue;
                }
            }

            if(currentMols.contains(p)){
                continue;
            }
            String cluster = oechem.OEGetSDData(p.getMol(),clusterTag);
            if(clustersList.contains(cluster)){
                continue;
            }
            return p;
        }
        return null;
    }

    public static Vector<PropertyMolecule> getDiverseMolecules(Vector<PropertyMolecule> molecules, int numCmpds){
        Random random = new Random(System.currentTimeMillis());
        Vector<PropertyMolecule> selectedMolecules = new Vector<PropertyMolecule>();
        int idx = random.nextInt(molecules.size());
        OEGraphMol firstMol = molecules.get(idx).getMol();
        OEFPDatabase fpdb = new OEFPDatabase(OEFPType.Path);
        fpdb.SetSimFunc(OESimMeasure.Tanimoto);
        fpdb.ClearCutoff();
        fpdb.AddFP(firstMol);
        Vector<Integer> selectedIdxes = new Vector<Integer>();
        selectedIdxes.add(idx);
        while(selectedMolecules.size()<numCmpds){
            double bestScore = 999999;
            int selectedIdx = -1;
            for(int i=0;i<molecules.size();i++){
                if(selectedIdxes.contains(i)){
                    continue;
                }
                PropertyMolecule propertyMolecule = molecules.get(i);
                double score1 = 0;
                OESimScoreIter oeSimScores = fpdb.GetScores(propertyMolecule.getMol());
                while(oeSimScores.hasNext()){
                    OESimScore simScore = oeSimScores.next();
                    score1 += simScore.GetScore();
                }
                if(score1<bestScore){
                    bestScore = score1;
                    selectedIdx = i;
                }
            }
            if(selectedIdx==-1){
                break;
            }
            PropertyMolecule bestMol = molecules.get(selectedIdx);
            selectedMolecules.add(bestMol);
            fpdb.AddFP(bestMol.getMol());
            selectedIdxes.add(selectedIdx);
        }
        return selectedMolecules;
    }

}
