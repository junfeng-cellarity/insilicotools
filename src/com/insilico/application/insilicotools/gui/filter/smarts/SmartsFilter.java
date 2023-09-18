package com.insilico.application.insilicotools.gui.filter.smarts;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import openeye.oechem.OEMatchBase;
import openeye.oechem.OESubSearch;
import openeye.oechem.oechem;

import java.util.ArrayList;

public class SmartsFilter extends Filter {
	private FilterResult filter(final PropertyMolecule[] molecules, final FilterState state) throws Exception {
		final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
		final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

		final SmartsFilterState filterState = (SmartsFilterState) state;
		OESubSearch subsearch = new OESubSearch();
		subsearch.Init(filterState.getSmarts());

		if (molecules.length > 0) {
			for (int i = 0; i < molecules.length; i++) {
				oechem.OEPrepareSearch(molecules[i].getMol(),subsearch);
				int matchCount = 0;
				for(OEMatchBase ma:subsearch.Match(molecules[i].getMol())){
					if(ma.IsValid()){
						matchCount += 1;
					}
				}
				if (matchCount>0) {
					if (!filterState.isMoreThanOne()) {
						pass.add(molecules[i]);
					} else if (matchCount > 1) {
						pass.add(molecules[i]);
					} else {
						fail.add(molecules[i]);
					}
				} else {
					fail.add(molecules[i]);
				}
			}
		}

		return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Matches", FilterSection.NEITHER), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]), "Does not match", FilterSection.NEITHER)});
	}

	public FilterResult filter(ProgressReporter reporter, PropertyMolecule[] mols, FilterState state) throws Exception {
		return filter(mols, state);
	}
}