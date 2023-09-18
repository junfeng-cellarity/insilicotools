package com.insilico.application.insilicotools.data;

/**
 * Created by jfeng1 on 4/27/17.
 */
public class CNSTempoObj extends MolProperty{
    private int chance;
    public CNSTempoObj(String property) {
        super(property);
        if(!isNumerical){
            isNumerical = true;
            value = 0;
            chance = 0;
            return;
        }
        if(value<1){
            chance = 3;
        }else if(value<2){
            chance = 8;
        }else if(value<3){
            chance = 17;
        }else if(value<4){
            chance = 39;
        }else if(value<5){
            chance = 53;
        }else{
            chance = 80;
        }

    }

    @Override
    public String toString() {
        return String.format("%s (%d%% CNS)",property, chance);
    }
}
