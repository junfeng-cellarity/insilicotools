package com.insilico.application.insilicotools.util;

import chemaxon.formats.MolExporter;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.deprotect.DeprotectionGroup;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;
import openeye.oechem.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.util.Vector;

/**
 * Created by jfeng1 on 2/1/16.
 */
public class Enumerator {
    OELibraryGen libraryGen;
    boolean hasReaction = false;
    String smirks;
    public Enumerator() {
    }

    public int getNumReactants(){
        return libraryGen.NumReactants();
    }

    public boolean isHasReaction() {
        return hasReaction;
    }

    public int setReagent(PropertyMolecule mol, int position){
        if(libraryGen!=null){
            OEGraphMol oemol = mol.getMol();
            oechem.OEAssignAromaticFlags(oemol);
            return libraryGen.SetStartingMaterial(oemol,position);
        }
        return 0;
    }

    public void clearReagents(){
        if(libraryGen!=null&&getNumReactants()>0) {
            for (int i = 0; i < getNumReactants(); i++) {
                libraryGen.ClearStartingMaterial(i);
            }
        }
    }

    public int addReagent(Vector<PropertyMolecule> mols, int position){
        if(libraryGen!=null){
            int matches = 0;
            for(PropertyMolecule mol:mols) {
                OEGraphMol oemol = mol.getMol();
                oechem.OEAssignAromaticFlags(oemol);
                matches += libraryGen.AddStartingMaterial(oemol, position);
            }
            return matches;
        }
        return 0;
    }

    public void enumerateToFile(oemolostream ofs, boolean deprotect, ProgressReporter progressReporter) throws IOException {
        OEMolBaseIter prodIter = libraryGen.GetProducts();
        int n = 0;
        OELibraryGen deprotector = null;
        File tempFile = null;
        oemolostream tempfs = null;
        if(deprotect){
            tempFile = File.createTempFile("enumeration",".oeb");
            tempfs = new oemolostream();
            tempfs.open(tempFile.getAbsolutePath());
            DeprotectionGroup boc = DeprotectionGroup.getDeprotectionGroups()[0];
            deprotector = new OELibraryGen();
            deprotector.Init(boc.getSmirks());
            deprotector.SetExplicitHydrogens(false);
            deprotector.SetValenceCorrection(true);
            deprotector.SetRemoveUnmappedFragments(true);
            if(deprotector.IsValid()){
                deprotect = true;
            }else{
                deprotect = false;
            }
        }
        Vector<String> product_smiles_list = new Vector<>();
        while(prodIter.hasNext()){
            OEMolBase mol = prodIter.next();
            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
            oechem.OEClearSDData(mol);
            String smiles = oechem.OEMolToSmiles(mol);
            if(product_smiles_list.contains(smiles)){
                continue;
            }else {
                product_smiles_list.add(smiles);
            }
            if(deprotect){
                oechem.OEGenerate2DCoordinates(mol);
                oechem.OEWriteMolecule(tempfs, mol);
                if(progressReporter!=null){
                    progressReporter.reportProgress("Enumerating ", n);
                }
                n+=1;
            }else {
                oechem.OEGenerate2DCoordinates(mol);
                oechem.OEWriteMolecule(ofs, mol);
                if(progressReporter!=null){
                    progressReporter.reportProgress("Enumerating ", n);
                }
                n+=1;
            }
        }
        product_smiles_list.clear();
        int max_products = n;
        n = 0;
        if(deprotect){
            tempfs.close();
            oemolistream ifs = new oemolistream();
            ifs.open(tempFile.getAbsolutePath());
            OEGraphMol product = new OEGraphMol();
            n = 1;
            while(oechem.OEReadMolecule(ifs,product)){
                deprotector.AddStartingMaterial(product,0);
                if(n%100==0||n>max_products){
                    OEMolBaseIter deprotectIter = deprotector.GetProducts();
                    while(deprotectIter.hasNext()){
                        OEMolBase mol = deprotectIter.next();
                        String smiles = oechem.OEMolToSmiles(mol);
                        if(product_smiles_list.contains(smiles)){
                            continue;
                        }
                        product_smiles_list.add(smiles);
                        oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
                        oechem.OEClearSDData(mol);
                        oechem.OEGenerate2DCoordinates(mol);
                        oechem.OEWriteMolecule(ofs, mol);
                        if(progressReporter!=null){
                            progressReporter.reportProgress("Deprotecting ", n);
                        }
                    }
                    deprotector.ClearStartingMaterial(0);
                }
                n += 1;
            }
            tempfs.close();
        }
    }


    public Vector<PropertyMolecule> enumerate(){
        Vector<String> dict = new Vector<String>();
        Vector<PropertyMolecule> products = new Vector<PropertyMolecule>();
        OEMolBaseIter prodIter = libraryGen.GetProducts();
        while(prodIter.hasNext()){
            OEMolBase mol = prodIter.next();
            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol);
            oechem.OEClearSDData(mol);
            String smiles = oechem.OEMolToSmiles(mol);
            if(dict.contains(smiles)){
                continue;
            }else{
                dict.add(smiles);
            }
            oechem.OEGenerate2DCoordinates(mol);
            if(mol.GetTitle().endsWith("_")){
                mol.SetTitle(mol.GetTitle().substring(0,mol.GetTitle().length()-1));
            }
            if(mol.GetTitle().startsWith("_")){
                mol.SetTitle(mol.GetTitle().substring(1));
            }
            products.add(new PropertyMolecule(new OEGraphMol(mol)));
        }
        return products;
    }

    public String getSmirks() {
        return smirks;
    }

    public void setSmirks(String smirks){
        this.smirks = smirks;
        libraryGen = new OELibraryGen();
        libraryGen.Init(smirks);
        libraryGen.SetExplicitHydrogens(false);
        libraryGen.SetValenceCorrection(true);
        libraryGen.SetRemoveUnmappedFragments(true);
        if(libraryGen.IsValid()){
            hasReaction = true;
        }else{
            hasReaction = false;
        }
    }

//    public void setRxn(String rxn){
//        libraryGen = new OELibraryGen();
//        oemolistream ifs = new oemolistream();
//        ifs.openstring(rxn);
//        OEQMol reaction = new OEQMol();
//        int opt = OEMDLQueryOpts.ReactionQuery | OEMDLQueryOpts.SuppressExplicitH;
//        oechem.OEReadMDLReactionQueryFile(ifs,reaction,opt);
//        ifs.close();
//        libraryGen.Init(reaction);
//        libraryGen.SetExplicitHydrogens(false);
//        libraryGen.SetValenceCorrection(true);
//        libraryGen.SetRemoveUnmappedFragments(true);
//        if(libraryGen.IsValid()){
//            hasReaction = true;
//        }else{
//            hasReaction = false;
//        }
//    }

//    private void calculateCoreMW(String smirks){
//        try {
//            Molecule mol = MolImporter.importMol(smirks, "smarts");
//            RxnMolecule rxnMol = RxnMolecule.getReaction(mol).clone();
//            if(rxnMol.getProductCount()==1){
//                Vector<MolAtom> atomTobeRemoved = new Vector<MolAtom>();
//                Molecule product = rxnMol.getProduct(0);
//                for(MolAtom atm:product.getAtomArray()){
//                    if(atm.getAtomMap()>0){
//                        atm.setAtno(MolAtom.ANY);
//                        atomTobeRemoved.add(atm);
//                    }
//                }
//                for(MolAtom atm:atomTobeRemoved){
//                    product.removeAtom(atm);
//                }
//                OEGraphMol oemol = ChemFunc.getMolFromMolString(MolExporter.exportToFormat(product, "mol"),OEFormat.MDL);
//                oechem.OETheFunctionFormerlyKnownAsStripSalts(oemol);
//                coreMW = oechem.OECalculateMolecularWeight(oemol);
////                ChemFunc.calculateOEPropertySingleMol(pmol);
//                //corePSA = pmol.getPSA();
////                coreMW = pmol.getMW();
//            }
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }

    public static void main(String[] args) {
        JFrame frame = new JFrame();
        final MSketchPane sketcher = MarvinFactory.getCompoundSketcher();

        JPanel p = new JPanel(new BorderLayout());
        JPanel btnPanel = new JPanel();
        final JButton enumerateBtn = new JButton("Enumerate");
        enumerateBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Molecule mol = sketcher.getMol();
                try {
                    String rxn = MolExporter.exportToFormat(mol, "smarts");
                    Enumerator enumerator = new Enumerator();
                    enumerator.setSmirks(rxn);
                    System.out.println(enumerator.isHasReaction());
                    System.out.println(enumerator.getNumReactants());
                    System.out.println(rxn);
                } catch (IOException e1) {
                    e1.printStackTrace();
                }


            }
        });
        btnPanel.add(enumerateBtn);
        p.add(sketcher,BorderLayout.CENTER);
        p.add(btnPanel,BorderLayout.SOUTH);

        frame.getContentPane().add(p);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(new Dimension(1024,768));
        frame.setVisible(true);

    }

}
