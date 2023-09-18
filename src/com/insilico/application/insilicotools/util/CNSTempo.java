package com.insilico.application.insilicotools.util;

/**
 * Created by jfeng1 on 4/24/17.
 */
public class CNSTempo {
    public static final int LOGP_TYPE = 0;
    public static final int NUM_ROT_BONDS = 1;
    public static final int NUM_ARO_RINGS = 2;
    public static final int NUM_CHAINS = 3;
    public static final int NUM_H_ACCEPTORS = 4;
    public static final int CARBON_HETERO_RATIO = 5;
    public static final int NUM_BASIC_AMINES = 6;
    public static final int NUM_NON_CONJUGATED_RING_CARBON = 7;
    float QL;
    float PL;
    float PU;
    float QU;
    float coeff;

    public CNSTempo(float QL, float PL, float PU, float QU, float coeff) {
        this.QL = QL;
        this.PL = PL;
        this.PU = PU;
        this.QU = QU;
        this.coeff = coeff;
    }

    public float getPenalty(float p){
        float PL_QL = PL-QL;
        float QU_PU = QU-PU;
        float penalty = -1;
        if(p>=PL&&p<=PU){
            penalty = 0;
        }else if(p<PL){
           if(PL_QL==0){
               penalty = 0;
           }else{
               penalty = (PL-p)/PL_QL;
           }
        }else{ //p>PU
            if(QU_PU==0){
                penalty = 0;
            }else{
                penalty = (p-PU)/QU_PU;
            }
        }
        return penalty;
    }

    public float getPenaltyCoeff(float p){
        return coeff*getPenalty(p);
    }

    public static CNSTempo getCNSTempoByType(int type){
        switch(type){
            case LOGP_TYPE:
                return new CNSTempo(-0.308f,2.608f,4.466f,5.6f,0.2568f);
            case NUM_ROT_BONDS:
                return new CNSTempo(0,1,4,8,1.0717f);
            case NUM_ARO_RINGS:
                return new CNSTempo(0,1,2,3,1.4692f);
            case NUM_CHAINS:
                return new CNSTempo(1,2,4,7,1.3882f);
            case NUM_H_ACCEPTORS:
                return new CNSTempo(1,2,3,6,0.7760f);
            case CARBON_HETERO_RATIO:
                return new CNSTempo(1.143f,2.143f,4.5f,10.5f,1.6542f);
            case NUM_BASIC_AMINES:
                return new CNSTempo(0,1,1,2,3.2205f);
            case NUM_NON_CONJUGATED_RING_CARBON:
                return new CNSTempo(0,0,3,10,0.0777f);
            default:
                return null;
        }
    }
}
