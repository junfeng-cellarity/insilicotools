package com.insilico.application.insilicotools.data;

/**
 * Created by jfeng1 on 4/27/17.
 */
public class CNSMPObj extends MolProperty{
    private int chance;
    public CNSMPObj(String property) {
        super(property);
        if(!isNumerical){
            isNumerical = true;
            value = 0;
            chance = 0;
            return;
        }
        if(value<1){
            chance = 0;
        }else if(value<2){
            chance = 8;
        }else if(value<3){
            chance = 15;
        }else if(value<4){
            chance = 24;
        }else if(value<5){
            chance = 31;
        }else{
            chance = 38;
        }

    }

    @Override
    public String toString() {
        return String.format("%s (%d%% CNS)",property, chance);
    }
}
