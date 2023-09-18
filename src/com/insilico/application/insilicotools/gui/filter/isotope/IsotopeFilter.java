package com.insilico.application.insilicotools.gui.filter.isotope;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.ArrayList;
import java.util.regex.Pattern;

public class IsotopeFilter extends Filter {
    private final Pattern isotopes = Pattern.compile("\\[[\\d].*");

    private FilterResult filter(final PropertyMolecule[] molecules) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

        for(int i = 0; i < molecules.length; i++) {
            PropertyMolecule m = molecules[i];
            if(!isotopes.matcher(m.getSmiles()).find()) {
                pass.add(m);
            }
            else {
                fail.add(m);
            }
        }
        return new FilterResult(pass.toArray(new PropertyMolecule[pass.size()]), fail.toArray(new PropertyMolecule[fail.size()]));
    }

    public FilterResult filter(ProgressReporter reporter, PropertyMolecule[] mols, FilterState state) throws Exception {
        return filter(mols);
    }
}