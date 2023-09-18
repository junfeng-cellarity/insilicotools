package com.insilico.application.insilicotools.gui.filter.SDField;

import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.Vector;

public class SDFieldFilterState extends FilterState {
    private SDFieldCondition[] conditions;
    private Vector<String> propertyKeys;

    public void setPropertyKeys(Vector<String> propertyKeys) {
        this.propertyKeys = propertyKeys;
    }

    public SDFieldCondition[] getConditions() {
        return conditions;
    }

    public Vector<String> getPropertyKeys() {
        return propertyKeys;
    }

    public void setConditions(final SDFieldCondition[] conditions) {
        this.conditions = conditions;
    }

    public String isValidState() {
        if(conditions.length == 0){
            return "At least one condition must be specified";
        }
        boolean valid = true;
        for (int i = 0; i < conditions.length; i++) {
            SDFieldCondition condition = conditions[i];
            if(condition.getOperator().equals(SDFieldCondition.BETWEEN)){
                if(condition.getText()==null||condition.getText().trim().length()==0){
                    valid = false;
                }else if(condition.getText2()==null||condition.getText2().trim().length()==0){
                    valid = false;
                }
            }else if(condition.getOperator().equals(SDFieldCondition.EQUALS)){
                if(condition.getText()==null||condition.getText().trim().length()==0){
                    valid = false;
                }
            }
            else{
                if(condition.getText()==null||condition.getText().trim().length()==0){
                    valid = false;
                }
            }
            if(!valid){
                break;
            }
        }
        if(valid){
            return null;
        }else{
            return "At least one condition is not valid.";
        }
    }

}
