package com.insilico.application.insilicotools.gui.filter.Neutralizer;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;

import java.util.ArrayList;

import com.insilico.application.insilicotools.gui.filter.Filter;
import com.insilico.application.insilicotools.gui.filter.FilterResult;
import com.insilico.application.insilicotools.gui.filter.FilterSection;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.OEAtomBase;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;


public class NeutralizerFilter extends Filter {

    private FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] molecules) throws Exception {
        final ArrayList<PropertyMolecule> pass = new ArrayList<PropertyMolecule>();
        final ArrayList<PropertyMolecule> fail = new ArrayList<PropertyMolecule>();
        if(molecules.length > 0) {
            for(int i = 0; i < molecules.length; i++) {
                reporter.reportProgress("",100*i/molecules.length);
                OEGraphMol mol = new OEGraphMol(molecules[i].getMol());
                if(neutralizeMol(mol)){
                    pass.add(new PropertyMolecule(mol));
                }else{
                    fail.add(molecules[i]);
                }
            }
        }
        return new FilterResult(new FilterSection[]{new FilterSection(pass.toArray(new PropertyMolecule[pass.size()]), "Neutralized"), new FilterSection(fail.toArray(new PropertyMolecule[fail.size()]),"Failed")});
    }

    private boolean neutralizeMol(OEGraphMol mol){
        boolean isCorrect = true;
        for(OEAtomBase atm:mol.GetAtoms()){
            if(atm.GetFormalCharge()==0){
                continue;
            }
            atm.SetFormalCharge(0);
            if(atm.GetAtomicNum()==7 && atm.GetExplicitValence()>3){
                isCorrect = false;
            }
        }
        oechem.OEAssignMDLHydrogens(mol);
        oechem.OEAssignAromaticFlags(mol);
        return isCorrect;
    }

    public static void main(String[] args) {
        OEGraphMol mol = new OEGraphMol();
        oechem.OESmilesToMol(mol,"CC[N+](C)C(=O)[O-]");
        OEGraphMol mol2 = new OEGraphMol(mol);
        System.out.println("Before:"+oechem.OEMolToSmiles(mol));
        boolean result = new NeutralizerFilter().neutralizeMol(mol);
        if(result) {
            System.out.println("After:" + oechem.OEMolToSmiles(mol));
        }else{
            System.out.println("After:" + oechem.OEMolToSmiles(mol2));
        }
    }

    public FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] mols, final FilterState state) throws Exception {
        return filter(reporter, mols);
    }
}
