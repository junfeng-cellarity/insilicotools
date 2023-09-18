package com.insilico.application.insilicotools.gui.filter.uniq;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.ArrayList;

public class UniqFilter extends Filter {
    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules, final FilterState state) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();
        ArrayList<String> smileList = new ArrayList<String>();
        if(molecules.length > 0) {
            for(PropertyMolecule mol:molecules){
                if(smileList.contains(mol.getSmiles())){
                    fail.add(mol);
                }else{
                    pass.add(mol);
                    smileList.add(mol.getSmiles());
                }
            }
        }

        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Unique"), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]), "Duplicates")});
        //return new FilterResult(new FilterSection[]{new FilterSection((PropertyMolecule[]) pass.toArray(new PropertyMolecule[pass.size()]), "Unique")});
    }
}
