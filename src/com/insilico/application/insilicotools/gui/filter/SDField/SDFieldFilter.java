package com.insilico.application.insilicotools.gui.filter.SDField;

import com.insilico.application.insilicotools.data.MolProperty;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.ArrayList;

public class SDFieldFilter extends Filter {
    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] mols, final FilterState state) throws Exception {
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        SDFieldFilterState filterState = (SDFieldFilterState) state;
        float lowerLimit = Float.NaN;
        float upperLimit = Float.NaN;
        float limit = Float.NaN;

        String pattern = null;
        for(int i=0;i<mols.length;i++){
            reporter.reportProgress("",100*i/mols.length);
            boolean passed = true;
            for(int j=0;j<filterState.getConditions().length;j++){
                SDFieldCondition condition = filterState.getConditions()[j];
                String key = condition.getProperty();
                MolProperty property = mols[i].getProperty(key);
                if(property==null){
                    passed = false;
                    break;
                }
                String propertyValue = property.getProperty();
                System.out.println(String.format("Key = %s and Property = %s",key,propertyValue));
                if(propertyValue==null){
                    passed = false;
                    break;
                }

                System.out.println(String.format("Operator: %s",condition.getOperator()));
                if(condition.getOperator().equals(SDFieldCondition.BETWEEN)){
                    float value = 0;
                    try {
                        value = Float.parseFloat(propertyValue);
                        lowerLimit = Float.parseFloat(condition.getText());
                        upperLimit = Float.parseFloat(condition.getText2());
                    } catch (NumberFormatException e1) {
                        passed = false;
                    }
                    System.out.println(String.format("%s %s",lowerLimit,upperLimit));
                    if(passed){
                        if(value < lowerLimit || value > upperLimit){
                            passed = false;
                        }
                    }
                }else if(condition.getOperator().equals(SDFieldCondition.EQUALS)){
                    boolean isNumeric = true;
                    float value = 0;
                    try {
                        value = Float.parseFloat(propertyValue);
                        limit = Float.parseFloat(condition.getText());
                    } catch (NumberFormatException e1) {
                        isNumeric = false;
                    }
                    if(isNumeric){
                        if(Math.abs(value-limit)>0.00001){
                            passed = false;
                        }
                    }else {
                        pattern = condition.getText().trim();
                        if (!propertyValue.trim().equalsIgnoreCase(pattern)) {
                            passed = false;
                        }
                    }
                }
                else{
                    float value = 0;
                    try {
                        value = Float.parseFloat(propertyValue);
                        limit = Float.parseFloat(condition.getText());
                    } catch (NumberFormatException e1) {
                        passed = false;
                    }
                    if(passed){
                        if(Math.abs(value-limit)<0.00001){
                            passed = false;
                        }else if(condition.getOperator().equals(SDFieldCondition.GREATER_THAN)){
                            if(value<limit){
                                passed = false;
                            }
                        }else{
                            if(value>limit){
                                passed = false;
                            }
                        }
                    }
                }

                if(!passed){
                    break;
                }
            }
            if(passed){
                pass.add(mols[i]);
            }else{
                fail.add(mols[i]);
            }
        }
        return new FilterResult(pass.toArray(new PropertyMolecule[pass.size()]), fail.toArray(new PropertyMolecule[fail.size()]));
    }

}
