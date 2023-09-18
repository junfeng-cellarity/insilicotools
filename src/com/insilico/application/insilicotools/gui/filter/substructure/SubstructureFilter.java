package com.insilico.application.insilicotools.gui.filter.substructure;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;

import javax.swing.*;
import java.util.ArrayList;

public class SubstructureFilter extends Filter {
    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] mols, final FilterState state) throws Exception {
        OESubSearch subSearch = new OESubSearch();
        SubstructureFilterState state1 = (SubstructureFilterState)state;
        subSearch.Init(state1.getQuery());
        ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

        for(int i = 0; i < mols.length; i++) {
            oechem.OEPrepareSearch(mols[i].getMol(),subSearch);
            final float counter = i;
            SwingUtilities.invokeLater(new Runnable() {
                public void run() {
                    reporter.reportProgress("Searching...", (int) ((counter / mols.length) * 100));
                }
            });

            try {
                if(subSearch.SingleMatch(mols[i].getMol())) {
                    pass.add(mols[i]);
                }
                else {
                    fail.add(mols[i]);
                }
            }
            catch(Exception err) {
                err.printStackTrace();
                fail.add(mols[i]);
            }
        }

        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Matches"), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]), "No Matches")});
    }

}