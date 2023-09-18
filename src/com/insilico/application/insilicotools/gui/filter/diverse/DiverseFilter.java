package com.insilico.application.insilicotools.gui.filter.diverse;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.util.ChemFunc;

import java.util.ArrayList;

public class DiverseFilter extends Filter {
	private FilterResult filter(final PropertyMolecule[] molecules, final FilterState state) throws Exception {
		final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
		final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();

		final DiverseFilterState filterState = (DiverseFilterState) state;
		if (molecules!=null&&molecules.length > 0&&filterState.getNumOfDiverseCompounds()>0) {
			ChemFunc.pickDiverseNMolecules(molecules,filterState.getNumOfDiverseCompounds());
			for (int i = 0; i < molecules.length; i++) {
				if(molecules[i].isSelected()){
					pass.add(molecules[i]);
					molecules[i].setIsSelected(false);
				}else{
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