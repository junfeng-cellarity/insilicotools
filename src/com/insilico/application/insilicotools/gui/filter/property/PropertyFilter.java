package com.insilico.application.insilicotools.gui.filter.property;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.util.ChemFunc;

import java.util.*;

public class PropertyFilter extends Filter {

	//acid pka and basic pka are not calculated from property tasks
	public final static String[] PROPERTIES_USED_FOR_FILTER = {
			"hydrogen-bond acceptors",
			"hydrogen-bond donors",
			"molecular weight",
			"sum of formal charges",
			"number of rings",
			"rotatable bonds",
			"XLogP",
			"2d PSA"
	};

	public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules, final FilterState state) throws Exception {
        ChemFunc.calculateOEProperty(new Vector<PropertyMolecule>(Arrays.asList(molecules)));
		final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();
		final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();

		PropertyFilterState filterState = (PropertyFilterState) state;

		final PropertyCondition[] conditions = filterState.getConditions();

		if (molecules.length > 0) {
            for(PropertyMolecule mol:molecules){
                boolean keep = true;
                for (int j = 0; j < conditions.length; j++) {
                    PropertyCondition condition = conditions[j];
                    double value = mol.getProperty(condition.getProperty()).getValue();

                    double compareValue = 0.0;
                    double compareValue2 = 0.0;
                    boolean use1 = true;
                    boolean use2 = true;

                    Object operator = condition.getOperator();
                    if (!operator.equals(PropertyCondition.BETWEEN)) {
                        use2 = false;
                    }


                    try {
                        compareValue = Double.parseDouble(condition.getText());
                    } catch (NumberFormatException err) {
                        use1 = false;
                    }

                    try {
                        compareValue2 = Double.parseDouble(condition.getText2());
                    } catch (NumberFormatException err) {
                        use2 = false;
                        if (operator.equals(PropertyCondition.BETWEEN)) {
                            continue;
                        }
                    }

                    if (!use1) continue;


                    if (!use2) {
                        if (operator.equals(PropertyCondition.EQUALS)) {
                            if(value != compareValue){
                                keep = false;
                                break;
                            }
                        } else {
                            if (operator.equals(PropertyCondition.GREATER_THAN)) {
                                if(value<=compareValue){
                                    keep = false;
                                    break;
                                }
                            } else {
                                if (operator.equals(PropertyCondition.LESS_THAN)) {
                                    if(value>=compareValue){
                                        keep = false;
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        if(value>compareValue2||value<compareValue){
                            keep = false;
                            break;
                        }
                    }
                }

                if (keep) {
                    pass.add(mol);
                } else {
                    fail.add(mol);
                }
            }
        }

		return new FilterResult(pass.toArray(new PropertyMolecule[pass.size()]), fail.toArray(new PropertyMolecule[fail.size()]));
	}
}