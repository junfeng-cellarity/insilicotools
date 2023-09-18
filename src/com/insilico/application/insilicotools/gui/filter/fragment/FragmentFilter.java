package com.insilico.application.insilicotools.gui.filter.fragment;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import openeye.oechem.oechem;

import java.util.ArrayList;

public class FragmentFilter extends Filter {
    private FilterResult filter(final PropertyMolecule[] molecules) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

        for(int i = 0; i < molecules.length; i++) {
            PropertyMolecule m = molecules[i];
            int[] parts = new int[m.getMol().GetMaxAtomIdx()];
            if(oechem.OEDetermineComponents(m.getMol(),parts)==1){
                pass.add(m);
            }
            else {
                fail.add(m);
            }
        }

        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Single Molecule"), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]), "Multiple Fragments")});
    }

    public FilterResult filter(ProgressReporter reporter, PropertyMolecule[] mols, FilterState state) throws Exception {
        return filter(mols);
    }
}