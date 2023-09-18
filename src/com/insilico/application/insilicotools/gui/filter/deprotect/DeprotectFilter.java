package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.util.Enumerator;
import openeye.oechem.OESDDataPair;
import openeye.oechem.oechem;

import java.util.ArrayList;
import java.util.Vector;


public class DeprotectFilter extends Filter {
	public DeprotectFilter() {
	}


	public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules, FilterState state) throws Exception {
		final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
		final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();
		DeprotectFilterState filterState = (DeprotectFilterState) state;
		Vector<Enumerator> enumerators = new Vector<Enumerator>();
		for(DeprotectionGroup group:filterState.getDeprotectionGroups()){
			Enumerator enumerator = new Enumerator();
			enumerator.setSmirks(group.getSmirks());
			enumerators.add(enumerator);
		}

		if (molecules.length > 0) {
			int n = 0;
			for(PropertyMolecule m:molecules){
				reporter.reportProgress("Deprotecting ...",100*n/molecules.length);
				PropertyMolecule product = null;
				for(Enumerator enumerator:enumerators){
					if(product==null){
						enumerator.setReagent(m,0);
					}else {
						enumerator.setReagent(product, 0);
					}
					Vector<PropertyMolecule> products = enumerator.enumerate();
					if(products.size()==1){
						product = products.get(0);
					}
				}
				if(product!=null){
					product.setName(m.getName());
					for(OESDDataPair pair:oechem.OEGetSDDataPairs(m.getMol())){
						oechem.OESetSDData(product.getMol(),pair.GetTag(),pair.GetValue());
					}
					pass.add(product);
				}else{
					fail.add(new PropertyMolecule(m.getMol()));
				}
				n++;
			}
		}
//        ChemFunc.calculateOEProperty(pass);
//        ChemFunc.calculateOEProperty(fail);
		return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Deprotected"), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]), "Unprotected")});
	}

}
