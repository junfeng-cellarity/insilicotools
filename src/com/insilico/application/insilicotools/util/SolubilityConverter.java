package com.insilico.application.insilicotools.util;

public class SolubilityConverter implements UnitConverter{
    double logS =0.0;

    public void set_value(double logS) {
        this.logS = logS;
    }

    public SolubilityConverter() {
    }

    @Override
    public double convert() {
        if(logS!=0.0){
            return Math.pow(10,logS)*1e6;
        }
        return 0;
    }
}
