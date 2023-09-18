package com.insilico.application.insilicotools.gui.filter.alert;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.*;

public class AlertFilter extends Filter {
    public static final int PAINS_FILTER = 0;
    public static final int TOX_FILTER = 1;
    public static final int LIBRARY_GUIDELINE_FILTER = 2;

    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules, final FilterState state) throws Exception {
        AlertFilterState filterState = (AlertFilterState) state;

        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final TreeMap<AlertRule,ArrayList<PropertyMolecule>> sections = new TreeMap<AlertRule,ArrayList<PropertyMolecule>>(new Comparator<AlertRule>() {
			public int compare(AlertRule o1, AlertRule o2) {
                return o1.getDescription().compareTo(o2.getDescription());
            }
        });

        Vector<AlertRule> selectedRules = filterState.getSelectedRules();
        for (int i = 0; i < molecules.length; i++) {
            Vector<AlertRule> hits = new Vector<AlertRule>();
            for(AlertRule rule:selectedRules){
                if(!rule.isPassed(molecules[i].getMol())){
                    hits.add(rule);
                }
            }
            reporter.reportProgress("",i * 100 / molecules.length);
            if (hits.isEmpty()) {
                pass.add(molecules[i]);
            } else {
                ArrayList<PropertyMolecule> members = sections.get(hits.get(0));
                if (members == null) {
                    members = new ArrayList<PropertyMolecule>();
                    sections.put(hits.get(0), members);
                }
                members.add(molecules[i]);
            }
        }

        FilterSection[] resultSections = new FilterSection[sections.size() + 1];
        resultSections[0] = new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Pass", FilterSection.PASS);
 
        Iterator<Map.Entry<AlertRule,ArrayList<PropertyMolecule>>> itr = sections.entrySet().iterator();
         for (int i = 1; itr.hasNext(); i++) {
        	 Map.Entry<AlertRule, ArrayList<PropertyMolecule>> entry = itr.next();
            resultSections[i] = new FilterSection(entry.getValue().toArray(new PropertyMolecule[0]), entry.getKey().getIndex() + " - " + entry.getKey().getDescription(), FilterSection.FAIL);
        }

        return new FilterResult(resultSections);
    }

}