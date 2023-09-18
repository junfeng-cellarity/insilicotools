package com.insilico.application.insilicotools.gui.filter.SDField;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;
import com.insilico.application.insilicotools.gui.filter.FilterController;

import java.util.HashSet;
import java.util.Vector;
import java.util.Iterator;

import openeye.oechem.OEGraphMol;
import openeye.oechem.OESDDataIter;
import openeye.oechem.OESDDataPair;
import openeye.oechem.oechem;

public class SDFieldFilterFactory extends FilterFactory {
    Vector<String> propertyKeys;

    public SDFieldFilterFactory(PropertyMolecule[] mols) {
        super("SD Fields");
        updateProperties(mols);
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new SDFieldFilterController(getName(), treeFilter, propertyKeys);
    }

    public void updateProperties(PropertyMolecule[] mols) {
        if(mols==null||mols.length==0){
            return;
        }
        int n = Math.min(10,mols.length);
        HashSet<String> propertyKeyHash = new HashSet<String>();
        for (int i = 0; i < n; i++) {
            PropertyMolecule mol = mols[i];
            for(String propertyName:mol.getPropertyNames()){
                propertyKeyHash.add(propertyName);
            }
            OEGraphMol oemol = mol.getMol();
            OESDDataIter iter = oechem.OEGetSDDataPairs(oemol);
            for(OESDDataPair dp:iter){
                propertyKeyHash.add(dp.GetTag());
            }
        }
        propertyKeys = new Vector<String>();
        for (Iterator<String> stringIterator = propertyKeyHash.iterator(); stringIterator.hasNext();) {
            String key = stringIterator.next();
            propertyKeys.add(key);
        }
    }

}
