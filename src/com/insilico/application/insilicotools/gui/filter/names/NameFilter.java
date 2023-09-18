package com.insilico.application.insilicotools.gui.filter.names;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashSet;

public class NameFilter extends Filter {

    private FilterResult filter(final PropertyMolecule[] molecules, final FilterState state) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

        final NameFilterState filterState = (NameFilterState) state;

        HashSet<String> names = new HashSet<String>();

        BufferedReader br = new BufferedReader(new StringReader(filterState.getNames()));
        String line;
        while((line = br.readLine()) != null) {
            if(line.trim().length() > 0) {
                names.add(line.trim());
            }
        }
        for(int i = 0; i < molecules.length; i++) {
            if(names.contains(molecules[i].getName().trim())) {
                pass.add(molecules[i]);
            }
            else {
                fail.add(molecules[i]);
            }
        }

        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Matches", FilterSection.NEITHER), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]), "Does not match", FilterSection.NEITHER)});
    }

    public FilterResult filter(ProgressReporter reporter, PropertyMolecule[] mols, FilterState state) throws Exception {
        return filter(mols, state);
    }

}