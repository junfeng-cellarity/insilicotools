package com.insilico.application.insilicotools.gui.filter.SaltStrip;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import java.util.ArrayList;

import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

public class SaltStripFilter extends Filter {

    private FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        if (molecules.length > 0) {
            for (int i = 0; i < molecules.length; i++) {
                reporter.reportProgress("",100*i/molecules.length);
                OEGraphMol oemol = new OEGraphMol(molecules[i].getMol());
                boolean result = oechem.OETheFunctionFormerlyKnownAsStripSalts(oemol);
                if(result){
                    PropertyMolecule e = new PropertyMolecule(oemol);
                    pass.add(e);
                }else {
                    pass.add(molecules[i]);
                }
            }
//            ChemFunc.calculateOEProperty(pass);
        }
        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Saltless")});
    }

    public FilterResult filter(ProgressReporter reporter, PropertyMolecule[] mols, FilterState state) throws Exception {
        return filter(reporter, mols);
    }

}
