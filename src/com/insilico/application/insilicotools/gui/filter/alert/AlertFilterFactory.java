package com.insilico.application.insilicotools.gui.filter.alert;


import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

import java.util.Vector;

public class AlertFilterFactory extends FilterFactory {

    public AlertFilterFactory() {
        super("Structure Alert");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new AlertFilterController(getName(), treeFilter);
    }

	protected static Vector<AlertRule> getAlertRules(int type) {
        if(type == AlertFilter.PAINS_FILTER){
            return AlertRule.generatePainsRules();
        }
        else{
            Vector<AlertRule> rules = new Vector<AlertRule>();
            for(AlertRule rule:AlertRule.LIBRARY_ALERT_RULES){
                rules.add(rule);
            }
            return rules;
        }
    }
}