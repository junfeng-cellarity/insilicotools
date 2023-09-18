package com.insilico.application.insilicotools.cmdline;
import chemaxon.calculations.pka.PKaTrainingResult;
import chemaxon.marvin.Trainer;
import chemaxon.calculations.pka.PKaTrainingUtil;
import chemaxon.marvin.calculations.pKaPlugin;

import java.io.File;
import java.io.IOException;

/**
 * Created by jfeng1 on 3/17/16.
 */
public class pKa_Training {
    public static void main(String[] args) {
//        try {
//            PKaTrainingResult pKaTrainingResult = PKaTrainingUtil.loadCorrectionData(new File("/Users/jfeng1/Databases/ADME/pka_training/irak1.pkadata"));
//            pKaPlugin plugin = new pKaPlugin();
//            plugin.setDoublePrecision(2);
//            plugin.setBasicpKaLowerLimit(-5);
//            plugin.setAcidicpKaUpperLimit(25);
//            plugin.setCorrectionData(pKaTrainingResult);
//
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//
//      cxtrain pka -i chemaxon_pka_train_2016_03_17 /Users/jfeng1/Databases/ChemAxon_pKa_Training/chemaxon_pKa_training_298_training.sdf

        try {
            Trainer.main(args);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
