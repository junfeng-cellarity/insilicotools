package com.insilico.application.insilicotools.gui.filter.protonate;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.util.Protonator;
import openeye.oechem.OEGraphMol;

import java.util.ArrayList;

public class ProtonatorFilter extends Filter {
    private FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();
        if(molecules.length > 0) {
            for(int i = 0; i < molecules.length; i++) {
                reporter.reportProgress("",100*i/molecules.length);
                OEGraphMol mol = new OEGraphMol(molecules[i].getMol());
                if(Protonator.getInstance().protonate(mol)){
                    pass.add(new PropertyMolecule(mol));
                }else{
                    fail.add(molecules[i]);
                }
            }
        }
        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Protonated"), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]),"Failed")});
    }


    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] mols, final FilterState state) throws Exception {
        return filter(reporter, mols);
    }

}
