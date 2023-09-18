package com.insilico.application.insilicotools.util;

import chemaxon.calculations.clean.Cleaner;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import chemaxon.util.MolHandler;
import com.insilico.application.insilicotools.InSilicoToolOptions;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import openeye.oechem.*;
import openeye.oegrid.OEGridFileType;
import openeye.oegrid.OEScalarGrid;
import openeye.oegrid.oegrid;
import openeye.oeomega.OEOmega;
import openeye.oeomega.OEOmegaOptions;
import openeye.oequacpac.OECharges;
import openeye.oequacpac.oequacpac;
import openeye.oeshape.*;
import openeye.oespicoli.OESurface;
import openeye.oespicoli.oespicoli;
import openeye.oeszybki.*;
import openeye.oezap.OEZap;
import org.apache.xmlrpc.XmlRpcException;
import org.junit.Assert;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class OEChemFunc {

	NumberFormat nf = new DecimalFormat("##.####");

	final static double rmstol = 2.0;

	public enum OEChargeTypes {
		MMFF94s, AM1BCC, B3LYP, B3LYP_OPTIMIZE
	}

	private static boolean initialized = false;
	static OEChemFunc _this;
	private static ExecutorService executor = Executors.newFixedThreadPool(2);

	public static OEChemFunc getInstance() {
		if (_this == null) {
			_this = new OEChemFunc();
		}
		return _this;
	}

	public OEChemFunc() {
		if (!initialized) {
			System.out.println("Initializing OEChem ...");
			initialize();
		}
	}

	public static void initialize() {
		try {
			OEChemWebLicenseInstaller.loadOELicenseFromWeb();

//			oechem.OEUseJavaHeap(false);
		} catch (IOException e) {
			e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
			JOptionPane.showMessageDialog(null, e.getMessage());
		}
		initialized = true;
	}

	public OEGraphMol getMol2D(OEGraphMol mol){
		Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(mol);
		Cleaner.clean(molecule,2);
		if(molecule.getExplicitHcount()>0){
			MolHandler mh = new MolHandler(molecule);
			mh.removeHydrogens();
		}
		OEGraphMol new_mol = OEChemFunc.getInstance().convertChemAxonMol(molecule);
		oechem.OECopySDData(new_mol,mol);
		return new_mol;
	}

	public OEGraphMol getMol3D(OEGraphMol mol2d){
		OEWallTimer timer = new OEWallTimer();
		timer.Start();
		if(mol2d==null){
			return null;
		}
		if(mol2d.GetDimension()==3){
			OEGraphMol mol =  new OEGraphMol(mol2d);
			oechem.OEAddExplicitHydrogens(mol);
			return mol;
		}
		Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(mol2d);
		Cleaner.clean(molecule,3);
		OEGraphMol mol = OEChemFunc.getInstance().convertChemAxonMol(molecule);
		oechem.OEAddExplicitHydrogens(mol);
		return mol;
	}

	public String fixMolFromPdb(String ligandPdbfile) {
		OEGraphMol oemol = new OEGraphMol();
		oemolistream ifs = new oemolistream();
		ifs.SetFormat(OEFormat.PDB);
		ifs.openstring(ligandPdbfile);
		if (ifs.IsValid()) {
			oechem.OEReadMolecule(ifs, oemol);
			ifs.close();
			if (oemol.NumAtoms() > 999) {
				return null;
			}
			oechem.OEAssignAromaticFlags(oemol);
			oemolostream molOutStream = new oemolostream();
			molOutStream.openstring();
			molOutStream.SetFormat(OEFormat.MDL);
			oechem.OEWriteMolecule(molOutStream, oemol);
			return molOutStream.GetString();
		} else {
			ifs.close();
			return null;
		}
	}

	public boolean inject3DCoordinatesFromPdb(OEGraphMol mol2D, String pdbString) {
		OEGraphMol mol3D = getOEMolFromPDBString(pdbString);
		OEQMol qmol = new OEQMol(mol2D);
		oechem.OESuppressHydrogens(qmol);
		OESubSearch subSearch = new OESubSearch(qmol, OEExprOpts.AtomicNumber, 0);
		OEMatchBaseIter iterator = subSearch.Match(mol3D, true);
		while (iterator.hasNext()) {
			OEMatchBase match = iterator.next();
			if (match.NumAtoms() == qmol.NumAtoms()) {
				OEMatchPairAtomIter atomIterator = match.GetAtoms();
				while (atomIterator.hasNext()) {
					OEMatchPairAtom atompair = atomIterator.next();
					int pattern_idx = atompair.getPattern().GetIdx();
					int target_idx = atompair.getTarget().GetIdx();
					OEDoubleArray coord = new OEDoubleArray(3);
					mol3D.GetCoords(mol3D.GetAtom(new OEHasAtomIdx(target_idx)), coord);
					mol2D.SetCoords(mol2D.GetAtom(new OEHasAtomIdx(pattern_idx)), coord);
				}
			}
			return true;
		}
		return false;
	}

	public Vector<Molecule> convertGlideResultSDF(String glideSDF, int numSolutions) {
		System.out.println(glideSDF);
		Vector<Molecule> mols = new Vector<Molecule>();

		OEGraphMol mol = new OEGraphMol();
		oemolistream molInStream = new oemolistream();
		molInStream.SetFormat(OEFormat.SDF);
		molInStream.openstring(glideSDF);
		while (molInStream.IsValid()) {
			oechem.OEReadMolecule(molInStream, mol);
			oechem.OESuppressHydrogens(mol);
			oemolostream molOutStream = new oemolostream();
			molOutStream.openstring();
			molOutStream.SetFormat(OEFormat.MDL);
			oechem.OEWriteMolecule(molOutStream, mol);

			Molecule newMol = null;
			try {
				newMol = MolImporter.importMol(molOutStream.GetString());
			} catch (MolFormatException e) {
				e.printStackTrace();
				return null;
			}
			for (OESDDataIter iter = oechem.OEGetSDDataPairs(mol); iter.hasNext(); ) {
				OESDDataPair dp = iter.next();
				if (dp.GetTag().equals("r_i_docking_score")) {
					newMol.setProperty("Docking Score", dp.GetValue());
				}
				if (dp.GetTag().equals("VdW")) {
					newMol.setProperty("VdW", dp.GetValue());
				}
				if (dp.GetTag().equals("Elec")) {
					newMol.setProperty("Elec", dp.GetValue());
				}
				if (dp.GetTag().equals("Strain")) {
					newMol.setProperty("Strain", dp.GetValue());
				}
				if (dp.GetTag().equals("_Total")) {
					newMol.setProperty("_Total", dp.GetValue());
				}
			}
			mols.add(newMol);
			if (mols.size() == numSolutions) {
				break;
			}
		}

		return mols;
	}


	public String canSmi(OEGraphMol mol, boolean isomeric, boolean kekule) {
		oechem.OEFindRingAtomsAndBonds(mol);
		oechem.OEAssignAromaticFlags(mol, OEAroModel.OpenEye);
		int smiflag = OESMILESFlag.Canonical;
		if (isomeric) {
			smiflag |= OESMILESFlag.ISOMERIC;
		}
		if (kekule) {
			for (OEBondBase bond : mol.GetBonds(new OEIsAromaticBond())) {
				bond.SetIntType(5);
			}
			oechem.OECanonicalOrderAtoms(mol);
			oechem.OECanonicalOrderBonds(mol);
			oechem.OEClearAromaticFlags(mol);
			oechem.OEKekulize(mol);
		}
		String smi = oechem.OECreateSmiString(mol, smiflag);
		return smi;
	}

	public double generateMultipleConformers(OEMol mol, int nConfs, double energyWindow) {
		OEOmega omega = new OEOmega();
		OEOmegaOptions opts = new OEOmegaOptions();
		opts.SetMaxConfs(nConfs);
		opts.SetStrictFrags(true);
		opts.SetFromCT(true);
		opts.SetEnergyWindow(energyWindow);
		opts.SetBuildForceField(OEForceFieldType.MMFF94S);
		opts.SetSearchForceField(OEForceFieldType.MMFF94S);
		omega.SetOptions(opts);
		omega.call(mol);
		double lowestEnergy = 9999;
		int count = 0;
		for (OEConfBaseIter iter = mol.GetConfs(); iter.hasNext(); ) {
			OEGraphMol oemol = new OEGraphMol(iter.next());
			double energy = getEnergyOE(oemol);
			if (count == 0) {
				lowestEnergy = energy;
				count += 1;
			} else if (energy < lowestEnergy) {
				lowestEnergy = energy;
			}
		}
		System.out.println("Lowest Energy:" + lowestEnergy);
		return lowestEnergy;
	}

	public double generateMultipleConformersMMFFs(OEMol mol, int nConfs, String forcefield, double rms, double ewindow) {
		OEOmega omega = new OEOmega();
		OEOmegaOptions opts = omega.GetOptions();
		opts.SetBuildForceField(OEForceFieldType.MMFF94S);
		opts.SetSearchForceField(OEForceFieldType.MMFF94S);
		opts.SetMaxConfs(nConfs);
		opts.SetEnergyWindow(ewindow);
		opts.SetRMSThreshold(rms);
		opts.SetEnumRing(true);
		opts.SetFromCT(true);
		opts.SetEnumNitrogen(1);
		opts.SetStrictFrags(true);
		omega.call(mol);
		System.out.println(omega.GetOptions().GetBuildForceField());
		System.out.println(omega.GetOptions().GetMaxConfs());
		System.out.println(omega.GetOptions().GetEnergyWindow());
		System.out.println(omega.GetOptions().GetRMSThreshold());
		double lowestEnergy = 9999;
		int count = 0;
		for (OEConfBaseIter iter = mol.GetConfs(); iter.hasNext(); ) {
			OEGraphMol oemol = new OEGraphMol(iter.next());
			double energy = getEnergyOE(oemol);
			if (count == 0) {
				lowestEnergy = energy;
				count += 1;
			} else if (energy < lowestEnergy) {
				lowestEnergy = energy;
			}
		}
		return lowestEnergy;
	}

	public void getRestrictiveConformers(OEMol mol, OEGraphMol template, int nConfs) {
		OEOmega omega = new OEOmega();
		omega.GetOptions().SetBuildForceField(OEForceFieldType.MMFF94S);
		omega.GetOptions().SetSearchForceField(OEForceFieldType.MMFF94S);
		omega.GetOptions().SetMaxConfs(nConfs);
		omega.GetOptions().SetFromCT(true);
//        omega.SetEnumNitrogen(1);
		OEQMol qmol = new OEQMol(template);
//        qmol.BuildExpressions(OEExprOpts.AtomicNumber, OEExprOpts.BondOrder);
		qmol.BuildExpressions(OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
		OESubSearch pat = new OESubSearch(qmol);
		//omega.SetFixMol(template,OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
		omega.SetFixQuery(pat);
		omega.SetFixRMS(0.5);
		omega.call(mol);
	}



	public double getEnergyOE(OEGraphMol oemol) {
		OESzybki sz = new OESzybki();
		sz.SetForceFieldType(openeye.oeszybki.OEForceFieldType.MMFF94S);
		sz.SetRunType(OERunType.SinglePoint);
		OESzybkiResults result = new OESzybkiResults();
		if (!sz.call(oemol, result)) {
			return 0.0;
		} else {
			return result.GetTotalEnergy();
		}
	}


	public Vector<SuperpositionSolution> mcsFlex2Flex(OEGraphMol template, OEGraphMol target, int nConfsTemplate, int nConfsTarget, ProgressReporter progressReporter) {
		OEMol refMol = new OEMol(template);
		progressReporter.reportProgress("Generating template conformations ...",0);
		generateMultipleConformers(refMol, nConfsTemplate, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		progressReporter.reportProgress("Generating target conformations ...",0);
		OEMol fitMol = new OEMol(target);
		double lowestEnergy = generateMultipleConformers(fitMol, nConfsTarget, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
		int totalConfs = refMol.NumConfs();
		int count = 0;
		Vector<Vector<Integer>> mcsMatchingList = null;
		for (OEConfBaseIter iter = refMol.GetConfs(); iter.hasNext(); ) {
			int percent = count++ * 100 / totalConfs;
			progressReporter.reportProgress(String.format("%d percent completed.", percent), percent);
			refMol.SetActive(iter.next());
			for (OEConfBaseIter iter2 = fitMol.GetConfs(); iter2.hasNext(); ) {
				fitMol.SetActive(iter2.next());
				OEGraphMol m1 = new OEGraphMol(refMol);
				OEGraphMol oeGraphMol = new OEGraphMol(fitMol);
				OEGraphMol m2 = new OEGraphMol(oeGraphMol);
				if (mcsMatchingList == null) {
					mcsMatchingList = getMCSOE(m1, oeGraphMol);
				}
				OEGraphMol m3 = mcsMatch(m1, m2, mcsMatchingList);
				if (m3 == null) {
					continue;
				}
				double score = Double.parseDouble(oechem.OEGetSDData(m3,"Shape/Feature Score"));
				double strainEnergy = getEnergyOE(oeGraphMol) - lowestEnergy;
				oechem.OESetSDData(m3,"deltaEnergy", nf.format(strainEnergy));
				SuperpositionSolution newSolution = new SuperpositionSolution(m1, m3, strainEnergy, score);
				addNewSolutions(solutions, score, newSolution);
			}
		}
		return solutions;
	}
	/*
    public Vector<SuperpositionSolution> mcsRigid2Flex(Molecule template, Molecule target, ProgressReporter progressReporter) {
        progressReporter.reportProgress(DesignProgressMonitor.INDETERMINATE, "MCS based alignment ...");
        Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
        try {
            Vector<Molecule> results = ModelingTasks.mcsAlignment(template.toFormat("mol"), target.toFormat("mol"));
            for (int i = 0; i < results.size(); i++) {
                Molecule molecule = results.elementAt(i);
                double strainEnergy = Double.parseDouble(molecule.getProperty("MMFF_CONF_Energy"));
                double comboScore = 0.5 * Double.parseDouble(molecule.getProperty("comboscore"));
                molecule.setProperty("deltaEnergy", nf.format(strainEnergy));
                molecule.setProperty("Shape/Feature Score", nf.format(comboScore));
                SuperpositionSolution solution = new SuperpositionSolution(template, molecule, strainEnergy, comboScore);
                solutions.add(solution);
            }
            Collections.reverse(solutions);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return solutions;
    }
    */

	public Vector<SuperpositionSolution> mcsRigid2Flex(OEGraphMol template, OEGraphMol target, int nConfs, ProgressReporter progressReporter) {
		progressReporter.reportProgress("Generating conformations ...",0);
		Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
		OEMol fitMol = new OEMol(target);
		double lowestEnergy = generateMultipleConformers(fitMol, nConfs, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);
		int totalConfs = fitMol.NumConfs();
		int count = 0;
		progressReporter.reportProgress("Searching for MCS ...",1);
		Vector<Vector<Integer>> mcsMatchingList = null;
		for (OEConfBaseIter iter = fitMol.GetConfs(); iter.hasNext(); ) {
			int percent = count++ * 100 / totalConfs;
			progressReporter.reportProgress(String.format("%d percent completed.", percent), percent);
			OEGraphMol targetMol = new OEGraphMol(iter.next());
			if (mcsMatchingList == null) {
				mcsMatchingList = getMCSOE(template, targetMol);
			}
			OEGraphMol resultMol = colorMatchByListOE(template, targetMol, mcsMatchingList, rmstol);
			if (resultMol == null) {
				continue;
			}
			double strainEnergy = getEnergyOE(targetMol) - lowestEnergy;
			oechem.OESetSDData(resultMol,"deltaEnergy", nf.format(new Double(strainEnergy)));
			double score = getComboScoreOE(template, resultMol);
			SuperpositionSolution newSolution = new SuperpositionSolution(new OEGraphMol(template), resultMol, strainEnergy, score);
			addNewSolutions(solutions, score, newSolution);
		}
		return solutions;
	}

	public void addNewSolutions(Vector<SuperpositionSolution> solutions, double score, SuperpositionSolution newSolution) {
		if (solutions.size() < InSilicoToolOptions.SUPERIMPOSE_NUM_RESULTS) {
			solutions.add(newSolution);
			Collections.sort(solutions);
		} else {
			if (score > solutions.get(0).getShapeScore()) {
				solutions.remove(0);
				solutions.add(newSolution);
				Collections.sort(solutions);
			}
		}
	}


	public OEGraphMol mcsMatch(OEGraphMol template, OEGraphMol target, Vector<Vector<Integer>> matchingList) {
		if (matchingList == null) {
			return null;
		}
		OEGraphMol bestMol = colorMatchByListOE(template, target, matchingList, rmstol);
		return bestMol;
	}

	public OEGraphMol pickingMatch(OEGraphMol template, OEGraphMol target, Vector<Integer> matchingList) {
		Vector<Vector<Integer>> allLists = new Vector<Vector<Integer>>();
		allLists.add(matchingList);
		OEGraphMol resultOEMol = colorMatchByListOE(template, target, allLists, rmstol);
		return resultOEMol;
	}

	public Vector<SuperpositionSolution> pickingFlex2Flex(OEGraphMol template, OEGraphMol target, int nConfsTemplate, int nConfsTarget, Vector<Integer> matchingList, ProgressReporter progressReporter) {
		Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
		OEMol refMol = new OEMol(template);
		progressReporter.reportProgress("Generating template conformations ...",0);
		generateMultipleConformers(refMol, nConfsTemplate, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		progressReporter.reportProgress("Generating target conformations ...",0);
		OEMol fitMol = new OEMol(target);
		double lowestEnergy = generateMultipleConformers(fitMol, nConfsTarget, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		int totalConfs = refMol.NumConfs();
		int count = 0;
		for (OEConfBaseIter iter = refMol.GetConfs(); iter.hasNext(); ) {
			int percent = count++ * 100 / totalConfs;
			progressReporter.reportProgress(String.format("%d percent completed.", percent), percent);
			refMol.SetActive(iter.next());
			OEGraphMol refGraphMol = new OEGraphMol(refMol);
			for (OEConfBaseIter iter2 = fitMol.GetConfs(); iter2.hasNext(); ) {
				fitMol.SetActive(iter2.next());
				OEGraphMol oeGraphMol = new OEGraphMol(fitMol);
				OEGraphMol m3 = pickingMatch(refGraphMol, oeGraphMol, matchingList);
				if (m3 == null) {
					continue;
				}
				double strainEnergy = getEnergyOE(oeGraphMol) - lowestEnergy;
				oechem.OESetSDData(m3,"deltaEnergy", nf.format(strainEnergy));
				double score = Double.parseDouble(oechem.OEGetSDData(m3,"Shape/Feature Score"));
				SuperpositionSolution newSolution = new SuperpositionSolution(refGraphMol, m3, strainEnergy, score);
				addNewSolutions(solutions, score, newSolution);
			}
		}
		return solutions;
	}

	public Vector<SuperpositionSolution> pickingRigid2Flex(Molecule template, Molecule target, int nConfs, Vector<Integer> matchingList, ProgressReporter progressReporter) {
//        stats.info(DesignStatUtils.createEZOverlayPicking(OEChemFunc.class));
		OEGraphMol refMol = convertChemAxonMol(template);
		OEGraphMol fitMol = convertChemAxonMol(target);
		return pickingRigid2Flex(refMol, fitMol, nConfs, matchingList, progressReporter);
	}

	public Vector<SuperpositionSolution> pickingRigid2Flex(OEGraphMol template, OEGraphMol target, int nConfs, Vector<Integer> matchingList, ProgressReporter progressReporter) {
		Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
		OEMol fitMol = new OEMol(target);
		progressReporter.reportProgress("Generating conformations ...",0);


		HashMap<Integer, Integer> indexMap = new HashMap<Integer, Integer>();
		for (OEAtomBaseIter iter = fitMol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			atom.SetMapIdx(atom.GetIdx() + 1);
		}

		double lowestEnergy = generateMultipleConformers(fitMol, nConfs, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		for (OEAtomBaseIter iter = fitMol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			if (atom.GetMapIdx() > 0) {
				indexMap.put(atom.GetMapIdx() - 1, atom.GetIdx());
			}
		}

		Vector<Integer> newMatchingList = new Vector<Integer>();
		for (int i = 0; i < matchingList.size(); i++) {
			if (i % 2 == 0) {
				newMatchingList.add(matchingList.elementAt(i));
			} else {
				Integer oldIndex = matchingList.elementAt(i);
				newMatchingList.add(indexMap.get(oldIndex));
			}
		}

		int totalConfs = fitMol.NumConfs();
		int count = 0;
		for (OEConfBaseIter iter = fitMol.GetConfs(); iter.hasNext(); ) {
			int percent = count++ * 100 / totalConfs;
			progressReporter.reportProgress(String.format("%d percent completed.", percent), percent);
			OEConfBase conf = iter.next();
			fitMol.SetActive(conf);

			OEGraphMol oeGraphMol = new OEGraphMol(fitMol);
			oechem.OESuppressHydrogens(oeGraphMol);
			OEGraphMol resultMol = pickingMatch(template, oeGraphMol, newMatchingList);
			if (resultMol == null) {
				continue;
			}
			double strainEnergy = getEnergyOE(oeGraphMol) - lowestEnergy;
			oechem.OESetSDData(resultMol,"deltaEnergy", nf.format(strainEnergy));
			double score = Double.parseDouble(oechem.OEGetSDData(resultMol,"Shape/Feature Score"));
			SuperpositionSolution newSolution = new SuperpositionSolution(template, resultMol, strainEnergy, score);
			addNewSolutions(solutions, score, newSolution);
		}
		return solutions;
	}

	public OEGraphMol matchByList(OEGraphMol ref, OEGraphMol fit, Vector<Vector<Integer>> matchingList){
		OEGraphMol bestMol = null;
		double[] rmat = new double[9];
		double[] trans = new double[3];
		double bestRmsd = 999.0;
		for (int j = 0; j < matchingList.size(); j++) {
			OEMatch xmatch = new OEMatch();
			Vector<Integer> list = matchingList.get(j);
			for (int k = 0; k < list.size(); k += 2) {
				int idx = list.elementAt(k);
				int idx_2 = list.elementAt(k + 1);
				OEAtomBase atom1 = ref.GetAtom(new OEHasAtomIdx(idx));
				OEAtomBase atom2 = fit.GetAtom(new OEHasAtomIdx(idx_2));
				xmatch.AddPair(atom1, atom2);
			}
			double rmsd = oechem.OERMSD(ref, fit, xmatch, true, rmat, trans);
			if (rmsd > bestRmsd) {
				continue;
			}
			if(bestMol!=null){
				bestMol.delete();
			}
			System.out.println(getStringFromOEMol(ref));
			bestMol = new OEGraphMol(fit);
			oechem.OERotate(bestMol, rmat);
			oechem.OETranslate(bestMol, trans);
			bestRmsd = rmsd;
		}
		oechem.OESetSDData(bestMol,"rmsd_to_reference",new DecimalFormat("#.##").format(bestRmsd));
		return bestMol;
	}

	public OEGraphMol colorMatchByListOE(OEGraphMol ref, OEGraphMol fit, Vector<Vector<Integer>> matchingList, double rmsdTolerance) {
		double bestTanimotoScore = 0.0;
		OEGraphMol bestMol = null;
		double[] rmat = new double[9];
		double[] trans = new double[3];

		NumberFormat nf = new DecimalFormat("##.###");
		for (int j = 0; j < matchingList.size(); j++) {
			OEMatch xmatch = new OEMatch();
			Vector<Integer> list = matchingList.get(j);
			for (int k = 0; k < list.size(); k += 2) {
				int idx = list.elementAt(k);
				int idx_2 = list.elementAt(k + 1);
				OEAtomBase atom1 = ref.GetAtom(new OEHasAtomIdx(idx));
				OEAtomBase atom2 = fit.GetAtom(new OEHasAtomIdx(idx_2));
				xmatch.AddPair(atom1, atom2);
			}
			double rmsd = oechem.OERMSD(ref, fit, xmatch, true, rmat, trans);

			if (rmsd > rmsdTolerance) {
				continue;
			}

			oechem.OERotate(fit, rmat);
			oechem.OETranslate(fit, trans);
			double tanimotoScore = getComboScoreOE(ref, fit);
			if (tanimotoScore > bestTanimotoScore) {
				bestMol = new OEGraphMol(fit);
				bestTanimotoScore = tanimotoScore;
			}
		}
		if (bestMol != null) {
			oechem.OESetSDData(bestMol, "Shape/Feature Score", nf.format(bestTanimotoScore));
		}
		return bestMol;
	}


	public OEGraphMol colorOverlapRigid2Rigid(OEGraphMol template, OEGraphMol target,boolean shapeOnly) {
		OEBestOverlay best = new OEBestOverlay();
		OEMol refMol = new OEMol(template);
		OEMol fitMol = new OEMol(target);
		best.SetRefMol(refMol);
		if(shapeOnly){
			best.ClearColorForceField();
		}else {
			best.SetColorForceField(OEColorFFType.ImplicitMillsDean);
			best.SetColorOptimize(true);
		}

		OEBestOverlayResultsIter resiter = best.Overlay(fitMol);
		OEBestOverlayScoreIter scoreiter = new OEBestOverlayScoreIter();
		oeshape.OESortOverlayScores(scoreiter, resiter, new OEHighestTanimotoCombo());
		OEBestOverlayScore score;
		OEGraphMol outMol = null;
		NumberFormat nf = new DecimalFormat("##.###");
		for (; scoreiter.IsValid(); scoreiter.Increment()) {
			score = scoreiter.Target();
			outMol = new OEGraphMol(fitMol.GetConf(new OEHasConfIdx(score.getFitconfidx())));
			score.Transform(outMol);
//            oechem.OESetSDData(outMol,"deltaEnergy",nf.format(getEnergyOE(outMol)-getEnergyOE(target)));
			oechem.OESetSDData(outMol, "Shape/Feature Score", nf.format(score.GetTanimotoCombo() / 2.0));
			//oechem.OESetSDData(outMol,"TanimotoCombo",nf.format(score.GetTanimotoCombo()));
			//oechem.OESetSDData(outMol,"ColorTanimoto",nf.format(score.GetColorTanimoto()));
			break;
		}
		return outMol;
	}

	public Vector<SuperpositionSolution> colorOverlapRigid2Flex(OEGraphMol template, OEGraphMol target, int nConfs, double rmsd, double energy_window, ProgressReporter progressReporter,int nSolutions, boolean shapeOnly) {
		Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
		OEBestOverlay best = new OEBestOverlay();
		OEMol refMol = new OEMol(template);
//		OEMol fitMol = new OEMol(target);
		if(progressReporter!=null) {
			progressReporter.reportProgress("Generating conformation ...", DesignProgressMonitor.INDETERMINATE);
		}
		Vector<PropertyMolecule> inputMols = new Vector<>();
		inputMols.add(new PropertyMolecule(target));
		Vector<OEMol> mcMols = ChemFunc.generateConformers(inputMols, energy_window, rmsd);
		if(mcMols.isEmpty()){
			return new Vector<>();
		}
		OEMol fitMol = mcMols.get(0);
//		double lowestEnergy = generateMultipleConformersMMFFs(fitMol, nConfs, "mmffs", rmsd, energy_window);
		if(progressReporter!=null) {
			progressReporter.reportProgress("Overlapping and scoring ...", DesignProgressMonitor.INDETERMINATE);
		}
		best.SetRefMol(refMol);
		if(shapeOnly){
			best.ClearColorForceField();
		}else {
			best.SetColorForceField(OEColorFFType.ImplicitMillsDean);
			best.SetColorOptimize(true);
		}
		OEBestOverlayResultsIter resiter = best.Overlay(fitMol);
		OEBestOverlayScoreIter scoreiter = new OEBestOverlayScoreIter();
		oeshape.OESortOverlayScores(scoreiter, resiter, new OEHighestComboScore());
		OEBestOverlayScore score;
		OEGraphMol outMol = null;
		NumberFormat nf = new DecimalFormat("##.###");
		for (; scoreiter.IsValid(); scoreiter.Increment()) {
			score = scoreiter.Target();
			outMol = new OEGraphMol(fitMol.GetConf(new OEHasConfIdx(score.getFitconfidx())));
			outMol.SetTitle(fitMol.GetTitle());
			oechem.OECopySDData(outMol,fitMol);
			score.Transform(outMol);
			double strainEnergy = Double.parseDouble(oechem.OEGetSDData(outMol,"deltaEnergy"));
			//double strainEnergy = getEnergyOE(outMol) ;
			//oechem.OESetSDData(outMol, "deltaEnergy", nf.format(strainEnergy));
			//oechem.OESetSDData(outMol,"TanimotoCombo",nf.format(score.GetTanimotoCombo()));
			double shapeScore = score.GetComboScore() / 2.0;
			oechem.OESetSDData(outMol, "Shape/Feature Score", nf.format(shapeScore));
			oechem.OESetSDData(outMol, "ColorTanimoto", nf.format(score.GetColorTanimoto()));
			SuperpositionSolution newSolution = new SuperpositionSolution(template, outMol, strainEnergy, shapeScore);
			solutions.insertElementAt(newSolution, 0);
			if (solutions.size() > 2*nSolutions) {
				break;
			}
		}
		Collections.sort(solutions, new Comparator<SuperpositionSolution>() {
			@Override
			public int compare(SuperpositionSolution o1, SuperpositionSolution o2) {
				Double score2 = o2.getScore();
				Double score1 = o1.getScore();
				if(Math.abs(score2-score1)<=0.01){
					return Double.compare(o1.strainEnergy, o2.strainEnergy);
				}
				return score2.compareTo(score1);
			}
		});
		Vector<SuperpositionSolution> results = new Vector<>();
		for(int i=0;i<Math.min(solutions.size(), nSolutions);i++){
			results.add(solutions.get(i));
		}
		return results;
	}

	public Vector<SuperpositionSolution> colorOverlapFlex2Flex(Molecule template, Molecule target, int nConfsTemplate, int nConfsTarget, ProgressReporter progressReporter) {
//        stats.info(DesignStatUtils.createEZOverlayShape(OEChemFunc.class));
		OEGraphMol templateOE = convertChemAxonMol(template);
		OEGraphMol targetOE = convertChemAxonMol(target);
		return colorOverlapFlex2Flex(templateOE, targetOE, nConfsTemplate, nConfsTarget, progressReporter);
	}

	public Vector<SuperpositionSolution> colorOverlapFlex2Flex(OEGraphMol template, OEGraphMol target, int nConfsTemplate, int nConfsTarget, ProgressReporter progressReporter) {
		Vector<SuperpositionSolution> solutions = new Vector<SuperpositionSolution>();
		OEMol mcMol = new OEMol(template);
		progressReporter.reportProgress("Generating template conformations ...",0);
		generateMultipleConformers(mcMol, nConfsTemplate, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		OEMol targetMol = new OEMol(target);
		progressReporter.reportProgress("Generating target conformations ...",0);
		double lowestEnergy = generateMultipleConformers(targetMol, nConfsTarget, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);

		OEBestOverlay best = new OEBestOverlay();
		best.SetColorForceField(OEColorFFType.ImplicitMillsDean);
		best.SetColorOptimize(true);

		best.SetRefMol(mcMol);
		OEBestOverlayResultsIter resiter = best.Overlay(targetMol);
		OEBestOverlayScoreIter scoreiter = new OEBestOverlayScoreIter();
		progressReporter.reportProgress("Sorting results ...",0);
		oeshape.OESortOverlayScores(scoreiter, resiter, new OEHighestComboScore());
		OEGraphMol resultMol = null;
		OEGraphMol templateMol = null;
		NumberFormat nf = new DecimalFormat("##.###");
		for (; scoreiter.IsValid(); scoreiter.Increment()) {
			OEBestOverlayScore score = scoreiter.Target();
			resultMol = new OEGraphMol(targetMol.GetConf(new OEHasConfIdx(score.getFitconfidx())));
			templateMol = new OEGraphMol(mcMol.GetConf(new OEHasConfIdx(score.getRefconfidx())));
			double strainEnergy = getEnergyOE(resultMol) - lowestEnergy;
			oechem.OESetSDData(resultMol, "deltaEnergy", nf.format(strainEnergy));
			//oechem.OESetSDData(resultMol,"TanimotoCombo",nf.format(score.GetTanimotoCombo()));
			double shapeScore = score.GetComboScore() / 2.0;
			oechem.OESetSDData(resultMol, "Shape/Feature Score", nf.format(shapeScore));
			oechem.OESetSDData(resultMol, "ColorTanimoto", nf.format(score.GetColorTanimoto()));
			score.Transform(resultMol);
			SuperpositionSolution newSolution = new SuperpositionSolution(templateMol, resultMol, strainEnergy, shapeScore);
			solutions.insertElementAt(newSolution, 0);
			if (solutions.size() == InSilicoToolOptions.SUPERIMPOSE_NUM_RESULTS) {
				break;
			}
		}
		return solutions;
	}


	public double getScaledColorScoreTanimotoOE(OEGraphMol oeTemplate, OEGraphMol oeTarget) {
		OEColorOverlap ov = new OEColorOverlap();
		ov.SetColorForceField(OEColorFFType.ImplicitMillsDean);
		ov.SetRefMol(oeTemplate);
		OEColorResults colorResults = new OEColorResults();
		ov.ColorScore(oeTarget, colorResults);
		return (double) colorResults.getScaledcolor();
	}


	public double getComboScoreOE(OEGraphMol template, OEGraphMol target) {
		double shapeScore = getAnalyticScoreTanimotoOE(template, target);
		double colorScore = getScaledColorScoreTanimotoOE(template, target);
		return (shapeScore + colorScore) / 2.0;
	}

	public double getAnalyticScoreTanimotoOE(OEGraphMol oeTemplate, OEGraphMol oeTarget) {
		OEOverlap ov = new OEOverlap();
		ov.SetMethod(OEOverlapMethod.Analytic);
		ov.SetRefMol(oeTemplate);

		OEOverlapResults res = new OEOverlapResults();
		ov.Overlap(oeTarget, res);
		return res.getTanimoto();
	}
/*
	public Molecule minimizeMolecule(Molecule molecule) {
		return minimizeMolecule(molecule, true);
	}

	public Molecule minimizeMolecule(Molecule molecule, boolean protonate) {
		OEGraphMol oemol = convertChemAxonMol(molecule);
		return convertOEChemMol(minimizeOEMol(oemol, null, protonate));
	}
*/

	public void removeHydrogenBrutal(Molecule molecule) {
		Vector<MolBond> bondToRemove = new Vector<MolBond>();
		Vector<MolAtom> atomToRemove = new Vector<MolAtom>();
		for (MolAtom atom : molecule.getAtomArray()) {
			for (int i = 0; i < atom.getBondCount(); i++) {
				MolBond bond = atom.getBond(i);
				MolAtom neighborAtom = bond.getOtherAtom(atom);
				if (neighborAtom.getAtno() == 1) {
					bondToRemove.add(bond);
					atomToRemove.add(neighborAtom);
				}
			}
		}
		for (MolBond bond : bondToRemove) {
			molecule.removeBond(bond);
		}

		for (MolAtom atom : atomToRemove) {
			molecule.removeAtom(atom);
		}
	}

	public Molecule fixHydrogen(Molecule molecule) {
		Molecule mol = molecule.cloneMolecule();
		removeHydrogenBrutal(mol);
		OEGraphMol oemol = convertChemAxonMol(mol);
		oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol);
		oechem.OESuppressHydrogens(oemol);
		oechem.OEAddExplicitHydrogens(oemol);
		Protonator.getInstance().protonate(oemol);
		Vector<Integer> atomIdList = new Vector<>();
		for(OEAtomBase atm:oemol.GetAtoms()){
			if(!atm.IsHydrogen()){
				atomIdList.add(atm.GetIdx());
			}
		}
		return convertOEChemMol(minimizeOEMol(oemol, atomIdList));
	}

	public OEGraphMol minimizeOEMol(OEGraphMol molecule, Vector<Integer> fixList){
		String molStr = getStringFromOEMol(molecule);
		if(fixList==null||fixList.isEmpty()){
			try {
				OEGraphMol mol = convertString2OEMol(ChemFunc.szybki_mol(molStr),OEFormat.SDF);
				oechem.OECopySDData(mol,molecule);
				return mol;
			} catch (MalformedURLException e) {
				e.printStackTrace();
			} catch (XmlRpcException e) {
				e.printStackTrace();
			}
		}else{
			StringBuilder sb = new StringBuilder();
			for(Integer id:fixList){
				sb.append(id);
				if(id!=fixList.get(fixList.size()-1)) {
					sb.append(",");
				}
			}
			try {
				OEGraphMol mol = convertString2OEMol(ChemFunc.szybki_mol_constrained(molStr,sb.toString()),OEFormat.SDF);
				return mol;
			} catch (MalformedURLException e) {
				e.printStackTrace();
			} catch (XmlRpcException e) {
				e.printStackTrace();
			}
		}
		return null;
	}

/*
	public OEGraphMol minimizeOEMol(OEGraphMol molecule, OEUnaryAtomPred pred, boolean protonate) {
		OESzybki sz = new OESzybki();
		sz.SetRunType(OERunType.CartesiansOpt);
		sz.SetForceFieldType(openeye.oeszybki.OEForceFieldType.MMFF94S);
		sz.SetOptCriteria(1000, 0.001);
		sz.SetOptimizerType(OEOptType.CG);

		OESzybkiResults result = new OESzybkiResults();

		if (pred != null) {
			sz.FixAtoms(pred);
		}
		if (protonate) {
			Protonator.getInstance().protonate(molecule);
		} else {
			Protonator.getInstance().adjust_hydrogens(molecule, false);
		}

		if (!sz.call(molecule, result)) {
			System.err.println("Failed to minimize structure.");
		} else {
			molecule.SetEnergy(result.GetTotalEnergy());
		}
		return molecule;
	}
 */

/*
	public OEGraphMol minimizeOEMol(OEGraphMol molecule, OEUnaryAtomPred pred, String pdbString) {
		if (pdbString == null) {
			return minimizeOEMol(molecule, pred, true);
		} else {
			OEGraphMol receptor = convertString2OEMol(pdbString, OEFormat.PDB);
			oechem.OESuppressHydrogens(receptor);
//          Protonator.getInstance().protonate(receptor);
			Protonator.getInstance().protonate(molecule);

			OESzybki sz = new OESzybki();
			sz.SetProtein(receptor);
			sz.SetRunType(OERunType.CartesiansOpt);
			sz.SetProteinDielectric(4.0);
			sz.SetProteinElectrostaticModel(OEProteinElectrostatics.ExactCoulomb);
			sz.SetForceFieldType(OEForceFieldType.MMFF94S);
			sz.SetOptimizerType(openeye.oeszybki.OEOptType.CG);
			OESzybkiResults result = new OESzybkiResults();
			if (pred != null) {
				sz.FixAtoms(pred);
			}

			if (!sz.call(molecule, result)) {
				return null;
			} else {
				molecule.SetEnergy(result.GetTotalEnergy() - result.GetInterEnergy());
				return molecule;
			}
		}
	}
*/

/*
	public Molecule minimizeWithinProtein(String ligandString, String pdbString) {
		//Load Receptor from String
		OEGraphMol receptor = convertString2OEMol(pdbString, OEFormat.PDB);

		//Load Ligand from String
		OEGraphMol ligand = convertString2OEMol(ligandString, OEFormat.MDL);
		Protonator.getInstance().protonate(ligand);

		OESzybkiResults res = new OESzybkiResults();

		OESzybki sz1 = new OESzybki();
		sz1.FixAtoms(new OEIsHeavy());
		sz1.SetForceFieldType(openeye.oeszybki.OEForceFieldType.MMFF94S);
		if (!sz1.call(ligand, res)) {
			System.err.println("failed to minimize ligand.");
			return null;
		}

		//double e1 = getEnergyOE(ligand);
//        double e1 = res.GetTotalEnergy();
//        System.out.println("Before minimization:" + e1);

		res.delete();
		res = new OESzybkiResults();
		OESzybki sz = new OESzybki();
		sz.SetProtein(receptor);
		sz.SetRunType(OERunType.CartesiansOpt);
		sz.SetProteinElectrostaticModel(OEProteinElectrostatics.ExactCoulomb);
		sz.SetProteinDielectric(4.0);
		sz.SetForceFieldType(openeye.oeszybki.OEForceFieldType.MMFF94S);
		sz.SetOptimizerType(openeye.oeszybki.OEOptType.CG);
		if (!sz.call(ligand, res)) {
			System.err.println("failed to minimize ligand.");
			return null;
		}

//        double vdw = res.GetEnergyTerm(OEPotentialTerms.VdWProteinLigand);
//        double es = res.GetEnergyTerm(OEPotentialTerms.CoulombProteinLigand);
		double e2 = res.GetTotalEnergy() - res.GetInterEnergy();

//        System.out.println("Vdw:"+vdw);
//        System.out.println("es"+es);
//        System.out.println("inter:"+res.GetInterEnergy());
//        System.out.println("total:"+res.GetTotalEnergy());
//        System.out.println("Strain"+e2);

		oemolostream molOutStream = new oemolostream();
		molOutStream.openstring();
		molOutStream.SetFormat(OEFormat.SDF);
//        oechem.OESuppressHydrogens(ligand);
		oechem.OEWriteMolecule(molOutStream, ligand);

		Molecule mol = null;
		try {
			mol = MolImporter.importMol(molOutStream.GetString());
		} catch (MolFormatException e) {
			return null;
		}
//        mol.setProperty("vdw",new Double(vdw).toString());
//        mol.setProperty("es",new Double(es).toString());
		mol.setPropertyObject("OEEnergy", new Double(e2));
		return mol;
	}
*/
	public Molecule getLowestEnergyStructure(Vector<Molecule> molVector) {
		if (molVector == null || molVector.isEmpty()) {
			return null;
		}
		if (molVector.size() == 1) {
			return molVector.get(0);
		}
		Molecule bestMol = molVector.get(0);
		double bestEnergy = Double.parseDouble(bestMol.getProperty("Energy"));
		for (int i = 1; i < molVector.size(); i++) {
			Molecule mol = molVector.elementAt(i);
			double energy = Double.parseDouble(mol.getProperty("Energy"));
			if (energy < bestEnergy) {
				bestMol = mol;
			}
		}
		return bestMol;
	}

	public Vector<PropertyMolecule> getMultiConformers(OEGraphMol mol, int nConfs, double ewindow, double rms, ProgressReporter pm) {
		Vector<PropertyMolecule> inputMols = new Vector<>();
		inputMols.add(new PropertyMolecule(mol));
		Vector<OEMol> mcMols = ChemFunc.generateConformers(inputMols, ewindow, rms);
		if (!mcMols.isEmpty()) {
			OEMol mcMol = mcMols.get(0);
			double lowestEnergy = 99999.0;
			int numConfs = 0;
			for (OEConfBaseIter iter = mcMol.GetConfs(); iter.hasNext(); iter.next()) {
				numConfs += 1;
			}

			if (pm != null) {
				pm.reportProgress("Generating conformers ...", -1);
			}

			Vector<OEGraphMol> oeMols = new Vector<OEGraphMol>();
			int count = 0;
			Vector<Future<String>> futures = new Vector<Future<String>>();
			for (OEConfBaseIter iter = mcMol.GetConfs(); iter.hasNext(); ) {
				OEConfBase conf = iter.next();
				OEGraphMol confMol = null;
				final String molString = ChemFunc.getMolString(new OEGraphMol(conf));
				if (molString != null && molString.length() > 0) {
					futures.add(executor.submit(new Callable<String>() {
						@Override
						public String call() throws Exception {
							//return ChemFunc.minimize_mol(molString);
							return molString; //no need to minimize after switch to macromodel
						}
					}));
				}
			}

			lowestEnergy = 99999;
			count = 0;
			for (Future<String> f : futures) {
				try {
					String minimized = f.get();
					if(pm!=null) {
						pm.reportProgress("Progress", 100 * (count++) / numConfs);
					}
					if (minimized == null) {
						continue;
					}

					OEGraphMol confMol = getOEMolFromStringWithData(minimized);
					if (confMol == null) {
						continue;
					}

					confMol.SetEnergy(Double.parseDouble(oechem.OEGetSDData(confMol, "r_mmod_Relative_Potential_Energy-S-OPLS"))/4.17828);

					if ((confMol.GetEnergy() - lowestEnergy) > ewindow) {
						continue;
					}

					if (confMol.GetEnergy() < lowestEnergy) {
						lowestEnergy = confMol.GetEnergy();
					}

					if (oeMols.isEmpty()) {
						oeMols.add(confMol);
					} else {
						OEGraphMol molToAdd = null;
						OEGraphMol molToRemove = null;
						boolean keep = true;
						for (OEGraphMol oemol : oeMols) {
							double rmsd = oechem.OERMSD(oemol, confMol, true, true, true);
							if (rmsd < rms) {
								if (confMol.GetEnergy() < oemol.GetEnergy()) {
									molToAdd = confMol;
									molToRemove = oemol;
								}
								keep = false;
								break;
							}
						}
						if (keep) {
							oeMols.add(confMol);
						}else{
							if(molToAdd!=null&&molToRemove!=null){
								oeMols.add(molToAdd);
								oeMols.remove(molToRemove);
							}
						}
					}
				} catch (Exception err) {
					err.printStackTrace();
				}
			}

			Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();

			Molecule tmpMol = convertOEChemMol(mol);
			MolHandler mh = new MolHandler(tmpMol);
			mh.removeHydrogens();
			OEGraphMol substructure = convertChemAxonMol(tmpMol);
			oechem.OESuppressHydrogens(substructure);
			count = 0;
			for(OEGraphMol oemol:oeMols){
				if (pm != null) {
					pm.reportProgress("Overlaying conformers...",100 * (count++) / numConfs);
				}
				oemol.SetEnergy(oemol.GetEnergy()-lowestEnergy);
				PropertyMolecule pmol = new PropertyMolecule(oemol);
				pmol.addProperty("Energy",nf.format(oemol.GetEnergy()));
				Vector<Vector<Integer>> subSearchMatchingList = OEChemFunc.getInstance().getSubSearchMatchingList(substructure,mol,oemol);
				if (subSearchMatchingList != null && subSearchMatchingList.size() > 0) {
					OEGraphMol fittedMol = OEChemFunc.getInstance().matchByList(mol, oemol, subSearchMatchingList);
					if (fittedMol != null) {
						pmol.setOEMol3D(fittedMol);
						pmol.addProperty("rmsd_to_reference",oechem.OEGetSDData(fittedMol,"rmsd_to_reference"));
					}
				}
				propertyMolecules.add(pmol);
			}
			Collections.sort(propertyMolecules, new Comparator<PropertyMolecule>() {
				public int compare(PropertyMolecule molecule, PropertyMolecule molecule1) {
					Double e1 = molecule.getProperty("Energy").getValue();
					Double e2 = molecule1.getProperty("Energy").getValue();
					return e1.compareTo(e2);
				}
			});

			if (propertyMolecules.size() > nConfs) {
				Vector<PropertyMolecule> propertyMolecules1 = new Vector<PropertyMolecule>(propertyMolecules.subList(0, nConfs - 1));
				propertyMolecules1.insertElementAt(new PropertyMolecule(mol), 0);

				return propertyMolecules1;
			} else {
				propertyMolecules.insertElementAt(new PropertyMolecule(mol), 0);
				return propertyMolecules;
			}

		}
		return new Vector<>();
	}

/*
	public Vector<PropertyMolecule> getMultiConformersOld(OEGraphMol mol, int nConfs, String forcefield, double rms, double ewindow, boolean protonate, ProgressReporter pm) {
		OEMol mcMol = new OEMol(mol);
		int numRotors = oechem.OECount(mol, new OEIsRotor());

        generateMultipleConformersMMFFs(mcMol, Math.max(numRotors * 50, nConfs), "mmffs", 0.1, 999.0);//generate as many as possible
        double lowestEnergy = 99999.0;
		int numConfs = 0;
		for (OEConfBaseIter iter = mcMol.GetConfs(); iter.hasNext(); iter.next()) {
			numConfs += 1;
		}

		if (pm != null) {
			pm.reportProgress("Generating conformers ...",0);
		}

		Vector<OEGraphMol> oeMols = new Vector<OEGraphMol>();
		int count = 0;
		Vector<Future<String>> futures = new Vector<Future<String>>();
		for (OEConfBaseIter iter = mcMol.GetConfs(); iter.hasNext(); ) {
			OEConfBase conf = iter.next();
			OEGraphMol confMol = null;
			if (pm != null) {
				if (forcefield.equals("mmffs")) {
					pm.reportProgress("Generating conformers...",100 * (count++) / numConfs);
				}
			}

			if (forcefield.equals("opls3")) {
				final String molString = getStringFromOEMol(new OEGraphMol(conf));
				if (molString != null && molString.length() > 0) {
					futures.add(executor.submit(new Callable<String>() {
						@Override
						public String call() throws Exception {
							return ChemFunc.minimize_mol(molString);
						}
					}));
				}
			} else {
				confMol = minimizeOEMol(new OEGraphMol(conf), null, protonate);
				if ((confMol.GetEnergy() - lowestEnergy) > ewindow) {
					continue;
				}
				if (confMol.GetEnergy() < lowestEnergy) {
					lowestEnergy = confMol.GetEnergy();
				}

				if (oeMols.isEmpty()) {
					oeMols.add(confMol);
				} else {
					boolean keep = true;
					OEGraphMol molToAdd = null;
					OEGraphMol molToRemove = null;
					for (OEGraphMol oemol : oeMols) {
						double rmsd = oechem.OERMSD(oemol, confMol, true, true, true);
						System.out.println(rmsd);
						if (rmsd < rms) {
							if (confMol.GetEnergy() < oemol.GetEnergy()) {
								molToAdd = confMol;
								molToRemove = oemol;
							}
							keep = false;
							break;
						}
					}
					if (keep) {
						oeMols.add(confMol);
					}else{
						if(molToAdd!=null&&molToRemove!=null){
							oeMols.add(molToAdd);
							oeMols.remove(molToRemove);
						}
					}
				}
			}

		}

		if (forcefield.equals("opls3")) {
			lowestEnergy = 99999;
			count = 0;
			for (Future<String> f : futures) {
				try {
					String minimized = f.get();
					if(pm!=null) {
						pm.reportProgress("Progress", 100 * (count++) / numConfs);
					}
					if (minimized == null) {
						continue;
					}

					OEGraphMol confMol = getOEMolFromStringWithData(minimized);
					if (confMol == null) {
						continue;
					}

					confMol.SetEnergy(Double.parseDouble(oechem.OEGetSDData(confMol, "r_ff_Potential_Energy-OPLS3")));

					if ((confMol.GetEnergy() - lowestEnergy) > ewindow) {
						continue;
					}

					if (confMol.GetEnergy() < lowestEnergy) {
						lowestEnergy = confMol.GetEnergy();
					}

					if (oeMols.isEmpty()) {
						oeMols.add(confMol);
					} else {
						OEGraphMol molToAdd = null;
						OEGraphMol molToRemove = null;
						boolean keep = true;
						for (OEGraphMol oemol : oeMols) {
							double rmsd = oechem.OERMSD(oemol, confMol, true, true, true);
							if (rmsd < rms) {
								if (confMol.GetEnergy() < oemol.GetEnergy()) {
									molToAdd = confMol;
									molToRemove = oemol;
								}
								keep = false;
								break;
							}
						}
						if (keep) {
							oeMols.add(confMol);
						}else{
							if(molToAdd!=null&&molToRemove!=null){
								oeMols.add(molToAdd);
								oeMols.remove(molToRemove);
							}
						}
					}
				} catch (Exception err) {
					err.printStackTrace();
				}
			}

		}
		Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();

		Molecule tmpMol = convertOEChemMol(mol);
		MolHandler mh = new MolHandler(tmpMol);
		mh.removeHydrogens();
		OEGraphMol substructure = convertChemAxonMol(tmpMol);
		oechem.OESuppressHydrogens(substructure);
		count = 0;
		for(OEGraphMol oemol:oeMols){
			if (pm != null) {
				pm.reportProgress("Overlaying conformers...",100 * (count++) / numConfs);
			}
			oemol.SetEnergy(oemol.GetEnergy()-lowestEnergy);
			PropertyMolecule pmol = new PropertyMolecule(oemol);
			pmol.addProperty("Energy",nf.format(oemol.GetEnergy()));
			Vector<Vector<Integer>> subSearchMatchingList = OEChemFunc.getInstance().getSubSearchMatchingList(substructure,mol,oemol);
			if (subSearchMatchingList != null && subSearchMatchingList.size() > 0) {
				OEGraphMol fittedMol = OEChemFunc.getInstance().matchByList(mol, oemol, subSearchMatchingList);
				if (fittedMol != null) {
					pmol.setOEMol3D(fittedMol);
					pmol.addProperty("rmsd_to_reference",oechem.OEGetSDData(fittedMol,"rmsd_to_reference"));
				}
			}
			propertyMolecules.add(pmol);
		}
		Collections.sort(propertyMolecules, new Comparator<PropertyMolecule>() {
			public int compare(PropertyMolecule molecule, PropertyMolecule molecule1) {
				Double e1 = molecule.getProperty("Energy").getValue();
				Double e2 = molecule1.getProperty("Energy").getValue();
				return e1.compareTo(e2);
			}
		});

		if (propertyMolecules.size() > nConfs) {
			Vector<PropertyMolecule> propertyMolecules1 = new Vector<PropertyMolecule>(propertyMolecules.subList(0, nConfs - 1));
			propertyMolecules1.insertElementAt(new PropertyMolecule(mol), 0);

			return propertyMolecules1;
		} else {
			propertyMolecules.insertElementAt(new PropertyMolecule(mol), 0);
			return propertyMolecules;
		}
	}
*/
/*
	private Vector<Molecule> getMultiConformers(Molecule molecule, int nConfs, boolean protonate) {
		Vector<Molecule> molVector = new Vector<Molecule>();
		OEGraphMol mol = convertChemAxonMol(molecule);
		OEMol mcMol = new OEMol(mol);
		generateMultipleConformersMMFFs(mcMol, nConfs, "mmffs", 1.0, InSilicoToolOptions.CONFGEN_ENERGY_WINDOW);
		for (OEConfBaseIter iter = mcMol.GetConfs(); iter.hasNext(); ) {
			OEConfBase conf = iter.next();
			OEGraphMol confMol = minimizeOEMol(new OEGraphMol(conf), null, protonate);
			Molecule confAXMol = convertOEChemMol(confMol);
			if (confAXMol == null) {
				continue;
			}
			confAXMol.setProperty("Energy", nf.format(confMol.GetEnergy()));
			molVector.add(confAXMol);
		}

		return molVector;
	}
*/
	public Vector<Vector<Integer>> getMCS(Molecule query, Molecule target) {
		OEGraphMol oeQmol = convertChemAxonMol(query);
		OEGraphMol oeTmol = convertChemAxonMol(target);
		return getMCSOE(oeQmol, oeTmol);
	}

	public Vector<Vector<Integer>> getSubSearchMatchingList(OEGraphMol substructure, OEGraphMol template, OEGraphMol fit){
		OESubSearch subsearch = new OESubSearch(substructure, OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
		oechem.OEPrepareSearch(template,subsearch);
		if (subsearch.GetMaxMatches() == 0) {
			return null;
		}
		OEMatchBaseIter iterator1 = subsearch.Match(template);

		oechem.OEPrepareSearch(fit,subsearch);
		if (subsearch.GetMaxMatches() == 0) {
			return null;
		}
		OEMatchBaseIter iterator2 = subsearch.Match(fit);

		Vector<Vector<Integer>> matchingList = new Vector<Vector<Integer>>();
		while (iterator1.hasNext()&&iterator2.hasNext()) {
			Vector<Integer> singleMatch = new Vector<Integer>();
			OEMatchBase oeMatchBase1 = iterator1.next();
			OEMatchPairAtomIter matchIterator1 = oeMatchBase1.GetAtoms();
			OEMatchBase oeMatchBase2 = iterator2.next();
			OEMatchPairAtomIter matchIterator2 = oeMatchBase2.GetAtoms();
			while (matchIterator1.hasNext()&&matchIterator2.hasNext()) {
				OEMatchPairAtom oeMatchPairAtom1 = matchIterator1.next();
				OEMatchPairAtom oeMatchPairAtom2 = matchIterator2.next();
				singleMatch.add(oeMatchPairAtom1.getTarget().GetIdx());
				singleMatch.add(oeMatchPairAtom2.getTarget().GetIdx());
//                System.out.println(String.format("%d %d",oeMatchPairAtom.getPattern().GetIdx(),oeMatchPairAtom.getTarget().GetIdx()));
			}
			matchingList.add(singleMatch);
		}
		return matchingList;
	}

	public Vector<Vector<Integer>> getMCSOE(OEGraphMol oeQmol, OEGraphMol oeTmol) {
		OEMCSSearch mcs = new OEMCSSearch(oeQmol, OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds);
		mcs.SetMinAtoms(4);
		OEMatchBaseIter iterator = mcs.Match(oeTmol);
		if (mcs.GetMaxMatches() == 0) {
			return null;
		}
		Vector<Vector<Integer>> matchList = new Vector<Vector<Integer>>();
		while (iterator.hasNext()) {
			Vector<Integer> singleMatch = new Vector<Integer>();
			OEMatchBase oeMatchBase = iterator.next();
			OEMatchPairAtomIter matchIterator = oeMatchBase.GetAtoms();
			while (matchIterator.hasNext()) {
				OEMatchPairAtom oeMatchPairAtom = matchIterator.next();
				singleMatch.add(oeMatchPairAtom.getPattern().GetIdx());
				singleMatch.add(oeMatchPairAtom.getTarget().GetIdx());
//                System.out.println(String.format("%d %d",oeMatchPairAtom.getPattern().GetIdx(),oeMatchPairAtom.getTarget().GetIdx()));
			}
			matchList.add(singleMatch);
		}
		return matchList;
	}

	public Molecule getTruncatedMolByMCS(Molecule patternMol, Molecule originalMol) {
		if (patternMol == null || originalMol == null) {
			return null;
		}
		/*
        if (dockingMol != null) {
            MolHandler mh = new MolHandler(dockingMol.cloneMolecule());
            mh.addHydrogens();

            HashMap<Integer, Integer> matchList2 = getFirstMCS(convertChemAxonMol(patternMol), convertChemAxonMol(mh.getMolecule()));
            if (matchList2 == null || matchList2.size() < patternMol.getAtomCount()) {
                return null;
            }
        }
        */
		HashMap<Integer, Integer> matchList = getFirstMCS_ConstrainedDocking(convertChemAxonMol(patternMol), convertChemAxonMol(originalMol));
		if (matchList == null || matchList.size() < patternMol.getAtomCount()) {
			return null;
		}

		Molecule targetMol = patternMol.cloneMolecule();
		targetMol.setDim(3);
		for (int i = 0; i < targetMol.getAtomCount(); i++) {
			if (matchList.get(new Integer(i)) == null) {
				return null;
			}
			int idx = matchList.get(new Integer(i));
			MolAtom atom1 = targetMol.getAtom(i);
			MolAtom atom2 = originalMol.getAtom(idx);
			if (atom2 == null) {
				return null;
			}
			atom1.setXYZ(atom2.getX(), atom2.getY(), atom2.getZ());
		}
		/*
        MolHandler mh = new MolHandler(targetMol);
        mh.addHydrogens();
        */
		return targetMol;
	}


	public OEGraphMol convertChemAxonMol(Molecule molecule) {
		if (molecule == null) {
			return null;
		}
		molecule.dearomatize();
		String molString = molecule.toFormat("mol");
		OEGraphMol mol = getOEMolFromString(molString);
		mol.SetTitle(molecule.getName());
		oechem.OESetSDData(mol,"Name",molecule.getName());
		return mol;
	}

	public OEGraphMol getOEMolFromPDBString(String pdbString) {
		OEGraphMol mol = new OEGraphMol();
		oemolistream molInStream = new oemolistream();
		molInStream.SetFormat(OEFormat.PDB);
		molInStream.openstring(pdbString);
		oechem.OEReadMolecule(molInStream, mol);
		molInStream.close();
		oechem.OEAssignAromaticFlags(mol);
		return mol;
	}

	public OEGraphMol getOEMolFromSmiles(String smiles) {
		OEGraphMol mol = new OEGraphMol();
		oemolistream molInStream = new oemolistream();
		molInStream.SetFormat(OEFormat.SMI);
		molInStream.openstring(smiles);
		oechem.OEReadMolecule(molInStream, mol);
		molInStream.close();
		oechem.OEAssignAromaticFlags(mol);
		return mol;
	}

	public OEGraphMol getOEMolFromStringWithData(String molString) {
		OEGraphMol mol = new OEGraphMol();
		oemolistream molInStream = new oemolistream();
		molInStream.SetFormat(OEFormat.SDF);
		molInStream.openstring(molString);
		oechem.OEReadMolecule(molInStream, mol);
		molInStream.close();
		oechem.OEAssignAromaticFlags(mol);
		return mol;
	}


	public OEGraphMol getOEMolFromString(String molString) {
		OEGraphMol mol = new OEGraphMol();
		oemolistream molInStream = new oemolistream();
		molInStream.SetFormat(OEFormat.MDL);
		molInStream.openstring(molString);
		oechem.OEReadMolecule(molInStream, mol);
		molInStream.close();
		oechem.OEAssignAromaticFlags(mol);
		return mol;
	}

	public String getStringFromOEMol(OEGraphMol mol) {
		return getStringFromOEMol(mol, OEFormat.MDL);
	}

	public String getStringFromOEMolWithData(OEGraphMol mol) {
		return getStringFromOEMol(mol, OEFormat.SDF);
	}


	public String getSmilesFromOEMol(OEGraphMol mol) {
		return getStringFromOEMol(mol, OEFormat.SMI);
	}

	public String getStringFromOEMol(OEGraphMol mol, int format) {
		oemolostream molOutStream = new oemolostream();
		molOutStream.openstring();
		molOutStream.SetFormat(format);
		oechem.OEWriteMolecule(molOutStream, mol);
		return molOutStream.GetString();
	}

	public Molecule convertOEChemMol(OEGraphMol mol) {
		return convertOEChemMol(mol, true);
	}

	public Molecule convertOEChemMol(OEGraphMol mol, boolean keepHydrogen) {
		if (mol == null) {
			return null;
		}
		OEGraphMol tmpMol = new OEGraphMol(mol);
		if (!keepHydrogen) {
			oechem.OESuppressHydrogens(tmpMol);
		}
		String molString = getStringFromOEMol(tmpMol);
		Molecule newMol = null;
		try {
			newMol = MolImporter.importMol(molString);
//            MolHandler mh = new MolHandler(tmpMol);
//            mh.removeHydrogens();
			//newMol = mh.getMolecule();
			for (OESDDataIter iter = oechem.OEGetSDDataPairs(tmpMol); iter.hasNext(); ) {
				OESDDataPair dp = iter.next();
				newMol.setProperty(dp.GetTag(), dp.GetValue());
			}
			newMol.setPropertyObject("OEEnergy", mol.GetEnergy());
		} catch (MolFormatException e) {
			System.err.println(e.getMessage());
		}
		return newMol;
	}

	public synchronized byte[] electrostatic_map(float ep_in, float ep_out, float grid_spacing, float buffer, OEGraphMol mol) throws Exception {
		OEZap oezap = new OEZap();
		oezap.SetInnerDielectric(ep_in);
		oezap.SetOuterDielectric(ep_out);
		oezap.SetGridSpacing(grid_spacing);
		oezap.SetBoundarySpacing(buffer);
		oezap.SetMolecule(mol);
		OEScalarGrid grid = new OEScalarGrid();
		System.err.println("calculating grid ...");
		if (oezap.CalcPotentialGrid(grid)) {
			System.err.println("writing ...");
			oeosstream stream = new oeosstream();
			oegrid.OEWriteGrid(stream, grid, OEGridFileType.Grasp);
			stream.close();
			return stream.getByteArray();
		} else {
			return null;
		}
	}

	public OEGraphMol loadPartialCharges(OEGraphMol mol, OEChargeTypes chargeType) throws Exception {
		oechem.OEAddExplicitHydrogens(mol);
		oechem.OEAssignBondiVdWRadii(mol);
		if (chargeType == OEChargeTypes.MMFF94s) {
			oechem.OEMMFFAtomTypes(mol);
			oechem.OEMMFF94PartialCharges(mol);
			return mol;
		} else if (chargeType == OEChargeTypes.AM1BCC) {
			boolean noHydrogen = false;
			boolean debug = false;
			oequacpac.OEAssignPartialCharges(mol, OECharges.AM1BCCSPt, noHydrogen, debug);
			return mol;
		}else {
			return null;
		}
	}

	public byte[] calcElectrostaticMap(String molString, int format, boolean isReceptor) throws Exception {
		OEGraphMol mol = convertString2OEMol(molString, format);
		if (isReceptor) {
			loadPartialCharges(mol, OEChargeTypes.MMFF94s);
			return electrostatic_map(2.0f, 80.0f, 1.0f, 2.0f, mol);
		} else {
			loadPartialCharges(mol, OEChargeTypes.AM1BCC);
			return electrostatic_map(2.0f, 80.0f, 0.5f, 2.0f, mol);
		}

	}

	public Vector calcElectrostaticMapWithMol(OEGraphMol mol, float ep_in, float ep_out, float grid_spacing, float buffer, OEChargeTypes type) throws Exception {
		Vector result = new Vector();
		if (type == null) {
			byte[] grid = electrostatic_map(ep_in, ep_out, grid_spacing, buffer, mol);
			result.add(grid);
			return result;
		} else {
			OEGraphMol newMol = loadPartialCharges(mol, type);
			byte[] grid = electrostatic_map(ep_in, ep_out, grid_spacing, buffer, newMol);
			result.add(grid);
			result.add(newMol);
			return result;
		}
	}

	public OEGraphMol convertString2OEMol(String molString, int format) {
		OEGraphMol mol = new OEGraphMol();
		oemolistream molInStream = new oemolistream();
		molInStream.SetFormat(format);
		molInStream.openstring(molString);
		oechem.OEReadMolecule(molInStream, mol);
		return mol;
	}

	public String rotateBond(String molString, int bgnTorsionIdx, int bgnIdx, int endIdx, int endTorsionIdx, double delta_theta) throws Exception {
		OEGraphMol oemol = getOEMolFromString(molString);
		if (!isBondRotor(bgnIdx, endIdx, oemol)) {
			throw new Exception("Not a rotable bond.");
		}
		double deltaDegree = delta_theta * Math.PI / 180.0;
		oechem.OESetTorsion(oemol, oemol.GetAtom(new OEHasAtomIdx(bgnTorsionIdx)), oemol.GetAtom(new OEHasAtomIdx(bgnIdx)), oemol.GetAtom(new OEHasAtomIdx(endIdx)), oemol.GetAtom(new OEHasAtomIdx(endTorsionIdx)), deltaDegree);
		return getStringFromOEMol(oemol);
	}

	public Vector<Molecule> torsionDrive(Molecule mol, int bgnTorsionIdx, int bgnIdx, int endIdx, int endTorsionIdx, double delta_theta) {

		OEGraphMol oemol = convertChemAxonMol(mol);
		if (!isBondRotor(bgnIdx, endIdx, oemol)) return null;

		Vector<Molecule> molVector = new Vector<Molecule>();
		molVector.add(mol.cloneMolecule());
		//setRotationPredicate(oemol,bgnIdx,endIdx);
		int n = (int) (360 / delta_theta);
		NumberFormat nf = new DecimalFormat("##.###");
		//oemolostream mos = new oemolostream("debug_all.sdf"); //for debug
		for (int i = 0; i < n; i++) {
			OEGraphMol oemolClone = new OEGraphMol(oemol);
			double deltaDegree = delta_theta * i * Math.PI / 180.0;
			oechem.OESetTorsion(oemolClone, oemolClone.GetAtom(new OEHasAtomIdx(bgnTorsionIdx)), oemolClone.GetAtom(new OEHasAtomIdx(bgnIdx)), oemolClone.GetAtom(new OEHasAtomIdx(endIdx)), oemolClone.GetAtom(new OEHasAtomIdx(endTorsionIdx)), deltaDegree);

			//OEGraphMol oemolClone = OEMolBondRotate(oemol,bgnIdx,endIdx,delta_theta*i*Math.PI/180.0);
//			OEUnaryAtomPred pred1 = new OEOrAtom(new OEHasAtomIdx(bgnIdx), new OEHasAtomIdx(endIdx));
//			pred1 = new OEOrAtom(pred1, new OEHasAtomIdx(bgnTorsionIdx));
//			pred1 = new OEOrAtom(pred1, new OEHasAtomIdx(endTorsionIdx));

			Vector<Integer> idList = new Vector<>();
			idList.add(bgnIdx);
			idList.add(bgnTorsionIdx);
			idList.add(endTorsionIdx);
			idList.add(endIdx);

			minimizeOEMol(oemolClone, idList);
			/*for debug
            oemolClone.SetTitle(String.format("%s_%s",oemol.GetTitle(),nf.format(delta_theta*i)));
            oechem.OESetSDData(oemolClone,"Energy",nf.format(oemolClone.GetEnergy()));
            oechem.OESetSDData(oemolClone,"Dihedral",nf.format(delta_theta*i));
            */
			Molecule minimizedMol = convertOEChemMol(oemolClone);
			minimizedMol.setProperty("Energy", nf.format(oemolClone.GetEnergy()));
			minimizedMol.setProperty("Dihedral", nf.format(delta_theta * i));
			molVector.add(minimizedMol);
			//  oechem.OEWriteMolecule(mos,oemolClone); //for debug
		}
		//mos.close();//for debug
		return molVector;
	}

	public boolean isBondRotor(Molecule mol, int bgnIdx, int endIdx) {
		OEGraphMol oemol = convertChemAxonMol(mol);
		return isBondRotor(bgnIdx, endIdx, oemol);
	}

	public boolean isBondRotor(int bgnIdx, int endIdx, OEGraphMol oemol) {
		OEAtomBase head = oemol.GetAtom(new OEHasAtomIdx(bgnIdx));
		OEAtomBase end = oemol.GetAtom(new OEHasAtomIdx(endIdx));
		OEBondBase bond = oemol.GetBond(head, end);
		if (bond == null) {
			System.err.println(String.format("no bond between %d and %d", bgnIdx, endIdx));
			return false;
		}
		/*
        if(!bond.IsRotor()){
            System.err.println("This bond is not a rotor.");
            return null;
        }
        */
		if (bond.GetOrder() == 2) {
			System.err.println("This bond is not a rotor.");
			return false;
		}

		if (bond.IsInRing()) {
			System.err.println("This bond is in a ring.");
			return false;
		}
		return true;
	}

	private static double getDistance(double[] xyz1, double[] xyz2) {
		double sumX = xyz1[0] - xyz2[0];
		double sumY = xyz1[1] - xyz2[1];
		double sumZ = xyz1[2] - xyz2[2];
		return Math.sqrt(sumX * sumX + sumY * sumY + sumZ * sumZ);
	}

	private void setRotationPredicate(OEGraphMol mol, int bgnIdx, int endIdx) {
		for (OEAtomBaseIter iter = mol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			atom.SetMapIdx(0);
		}
		setRotationMark(mol, bgnIdx, endIdx, bgnIdx, endIdx);

		for (OEAtomBaseIter iter = mol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			System.out.println(String.format("%d : %d", atom.GetIdx(), atom.GetMapIdx()));
		}

	}

	private void setRotationMark(OEGraphMol mol, int bgnIdx, int endIdx, int origBgnIdx, int origEndIdx) {
		OEAtomBase bgnAtom = mol.GetAtom(new OEHasAtomIdx(bgnIdx));
		OEAtomBase endAtom = mol.GetAtom(new OEHasAtomIdx(endIdx));
		for (OEBondBaseIter iter = endAtom.GetBonds(); iter.hasNext(); ) {
			OEBondBase bond = iter.next();
			OEAtomBase nbr = bond.GetNbr(endAtom);
			if (nbr == bgnAtom) {
				continue;
			}
			if (nbr.GetIdx() == origBgnIdx || nbr.GetIdx() == origEndIdx) {
				continue;
			}
			if (nbr.GetMapIdx() == 0) {
				nbr.SetMapIdx(1);
				setRotationMark(mol, endAtom.GetIdx(), nbr.GetIdx(), origBgnIdx, origEndIdx);
			}
		}
	}

	public OEGraphMol OEMolBondRotate(OEGraphMol mol, int bgnIdx, int endIdx, double theta) {
		OEAtomBase head = mol.GetAtom(new OEHasAtomIdx(bgnIdx));
		OEAtomBase end = mol.GetAtom(new OEHasAtomIdx(endIdx));
		OEBondBase bond = mol.GetBond(head, end);
		if (bond == null || !bond.IsRotor() || bond.IsInRing()) {
			System.err.println(String.format("no bond between %d and %d", bgnIdx, endIdx));
			return null;
		}
		double[] headXyz = new double[3];
		double[] endXyz = new double[3];
		mol.GetCoords(head, headXyz);
		mol.GetCoords(end, endXyz);

		double headX = headXyz[0];
		double headY = headXyz[1];
		double headZ = headXyz[2];

		double endX = endXyz[0];
		double endY = endXyz[1];
		double endZ = endXyz[2];

		double bondLength = getDistance(headXyz, endXyz);
		OEGraphMol molClone = new OEGraphMol(mol);

		double cos_alpha_1 = (headX - endX) / bondLength;
		double cos_beta_1 = (headY - endY) / bondLength;
		double cos_gamma_1 = (headZ - endZ) / bondLength;
		for (OEAtomBaseIter iter = molClone.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			if (atom.GetMapIdx() == 0) {
				continue;
			}
			if (atom == head || atom == end) {
				continue;
			}
			double[] atomXyz = new double[3];
			molClone.GetCoords(atom, atomXyz);
			double atomX = atomXyz[0];
			double atomY = atomXyz[1];
			double atomZ = atomXyz[2];
			double d_ac = getDistance(atomXyz, endXyz);

			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			double cos_alpha = (atomX - endX) / d_ac;
			double cos_beta = (atomY - endY) / d_ac;
			double cos_gamma = (atomZ - endZ) / d_ac;
			double cos_sum = cos_alpha * cos_alpha_1 + cos_beta * cos_beta_1 + cos_gamma * cos_gamma_1;
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//                     if(cos_sum >= 1.0) cos_sum = 0.9999;
//                     if(cos_sum <= -1.0) cos_sum = -0.9999;
//                     assert(cos_sum < 1.0 && cos_sum > -1.0);

			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			double phi = Math.acos(cos_sum);
			double cos_alpha_2 = atomX - endX - d_ac * Math.cos(phi) * cos_alpha_1;
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

			cos_alpha_2 = cos_alpha_2 / (d_ac * Math.sin(phi));

			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			double cos_beta_2 = atomY - endY - d_ac * Math.cos(phi) * cos_beta_1;
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

			cos_beta_2 = cos_beta_2 / (d_ac * Math.sin(phi));

			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			double cos_gamma_2 = atomZ - endZ - d_ac * Math.cos(phi) * cos_gamma_1;
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

			cos_gamma_2 = cos_gamma_2 / (d_ac * Math.sin(phi));

			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
			double cos_alpha_3 = cos_beta_1 * cos_gamma_2 - cos_gamma_1 * cos_beta_2;
			double cos_beta_3 = cos_gamma_1 * cos_alpha_2 - cos_alpha_1 * cos_gamma_2;
			double cos_gamma_3 = cos_alpha_1 * cos_beta_2 - cos_beta_1 * cos_alpha_2;
			/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

			double newAtomX = endX + d_ac * Math.cos(phi) * cos_alpha_1 + d_ac * Math.sin(phi) * Math.cos(theta) * cos_alpha_2 + d_ac * Math.sin(phi) * Math.sin(theta) * cos_alpha_3;
			double newAtomY = endY + d_ac * Math.cos(phi) * cos_beta_1 + d_ac * Math.sin(phi) * Math.cos(theta) * cos_beta_2 + d_ac * Math.sin(phi) * Math.sin(theta) * cos_beta_3;
			double newAtomZ = endZ + d_ac * Math.cos(phi) * cos_gamma_1 + d_ac * Math.sin(phi) * Math.cos(theta) * cos_gamma_2 + d_ac * Math.sin(phi) * Math.sin(theta) * cos_gamma_3;

			double[] newXyz = new double[3];
			newXyz[0] = newAtomX;
			newXyz[1] = newAtomY;
			newXyz[2] = newAtomZ;

			molClone.SetCoords(atom, newXyz);
		}
		return molClone;
	}

	public HashMap<Integer, Integer> getFirstMCS(OEGraphMol oeQmol, OEGraphMol oeTmol, int type) {
		OEMCSSearch mcs = new OEMCSSearch(oeQmol, OEExprOpts.DefaultAtoms, OEExprOpts.DefaultBonds, type);
		mcs.SetMCSFunc(new OEMCSMaxBondsCompleteCycles());
		mcs.SetMinAtoms(4);
		OEMatchBaseIter iterator = mcs.Match(oeTmol);
		if (mcs.GetMaxMatches() == 0) {
			return null;
		}
		HashMap<Integer, Integer> singleMatch = null;
		if (iterator.hasNext()) {
			singleMatch = new HashMap<Integer, Integer>();
			OEMatchBase oeMatchBase = iterator.next();
			OEMatchPairAtomIter matchIterator = oeMatchBase.GetAtoms();
			while (matchIterator.hasNext()) {
				OEMatchPairAtom oeMatchPairAtom = matchIterator.next();
				singleMatch.put(oeMatchPairAtom.getPattern().GetIdx(), oeMatchPairAtom.getTarget().GetIdx());
			}
		}
		return singleMatch;
	}

	public HashMap<Integer, Integer> getFirstMCS_ConstrainedDocking(OEGraphMol oeQmol, OEGraphMol oeTmol) {
		OEMCSSearch mcs = new OEMCSSearch(oeQmol, OEExprOpts.AtomicNumber, OEExprOpts.BondOrder, OEMCSType.Exhaustive);
		mcs.SetMCSFunc(new OEMCSMaxAtomsCompleteCycles());
		mcs.SetMinAtoms(4);
		OEMatchBaseIter iterator = mcs.Match(oeTmol);
		if (mcs.GetMaxMatches() == 0) {
			return null;
		}
		HashMap<Integer, Integer> singleMatch = null;
		if (iterator.hasNext()) {
			singleMatch = new HashMap<Integer, Integer>();
			OEMatchBase oeMatchBase = iterator.next();
			OEMatchPairAtomIter matchIterator = oeMatchBase.GetAtoms();
			while (matchIterator.hasNext()) {
				OEMatchPairAtom oeMatchPairAtom = matchIterator.next();
				singleMatch.put(oeMatchPairAtom.getPattern().GetIdx(), oeMatchPairAtom.getTarget().GetIdx());
			}
		}
		return singleMatch;
	}

	public static Vector<OEGraphMol> getChiralMols(OEGraphMol originalMol, boolean skipDefined){
		oechem.OEPerceiveChiral(originalMol);
		OEGraphMol mol = new OEGraphMol(originalMol);
		Vector<OEGraphMol> molList = new Vector<OEGraphMol>();
		Vector<Integer> chiralAtmIds = new Vector<Integer>();
		int numChiralCenters = 0;
		HashMap<Integer,Integer> chiral_id_dict = new HashMap<>();
		HashMap<Integer,String> chiral_dict = new HashMap<>();
		for (OEAtomBaseIter iter0 = mol.GetAtoms(); iter0.hasNext(); ) {
			OEAtomBase atm  = iter0.next();
			if(atm.IsChiral()&&atm.GetAtomicNum()==OEElemNo.C){
				numChiralCenters += 1;
				Integer atom_id = new Integer(atm.GetIdx());
				chiral_id_dict.put(numChiralCenters,atom_id);
				if(atm.HasStereoSpecified(OEAtomStereo.Tetrahedral)){
					boolean found = false;
					OEGroupBaseIter oeGroupBases = mol.GetGroups(new OEHasGroupType(OEGroupType.MDLAndStereo));
					for(OEGroupBase group:oeGroupBases){
						if(group.HasAtom(atm)){
							chiral_dict.put(atom_id,"rac");
							chiralAtmIds.add(atom_id);
							found = true;
							break;
						}
					}
					if(!found){
						chiral_dict.put(atom_id, "abs");
						if(!skipDefined) {
							chiralAtmIds.add(atom_id);
						}
					}
				}else{
					chiral_dict.put(atom_id,"rac");
                    chiralAtmIds.add(atom_id);
				}
			}
		}
		if(chiralAtmIds.isEmpty()){
			if(numChiralCenters>0){
				oechem.OESetSDData(mol,"No_Stereocenters",""+numChiralCenters);
				for(int i=0;i<numChiralCenters;i++) {
					int chiral_id = i+1;
					System.out.println(chiral_id);
					oechem.OESetSDData(mol, String.format("Stereocenter_%d",chiral_id),chiral_dict.get(chiral_id_dict.get(chiral_id)));
				}
			}
			molList.add(mol);
			return molList;
		}
		String[] rs = new String[]{"R","L"};
		String[] rsList = getAllLists(rs,chiralAtmIds.size());
		for(String rsConfig:rsList){
			OEGraphMol newMol = getChiralMol(mol,rsConfig,chiralAtmIds);
			if(numChiralCenters>0){
				oechem.OESetSDData(newMol,"No_Stereocenters",""+numChiralCenters);
				Vector<String> chirals = new Vector<>();
				for(int i=0;i<numChiralCenters;i++) {
					int chiral_id = i+1;
					String chiral = chiral_dict.get(chiral_id_dict.get(chiral_id));
					chirals.add(chiral);
//					oechem.OESetSDData(newMol, String.format("Stereocenter_%d",chiral_id), chiral);
				}
				Collections.sort(chirals, new Comparator<String>() {
					@Override
					public int compare(String o1, String o2) {
						if(o1.equals(o2)){
							return 0;
						}
						if(o1.equals("rac")){
							return -1;
						}
						if(o1.equals("abs")){
							return 1;
						}
						return 0;
					}
				});
				for(int i=0;i<chirals.size();i++){
					oechem.OESetSDData(newMol, String.format("Stereocenter_%d",i+1), chirals.get(i));
				}
			}
			molList.add(newMol);
		}
		return  molList;
	}

	private static OEGraphMol getChiralMol(OEGraphMol mol, String rsConfig, Vector<Integer> atomIds){
		OEGraphMol newMol = new OEGraphMol(mol);
		oechem.OEPerceiveChiral(newMol);
		Assert.assertTrue(rsConfig!=null&&atomIds!=null&&rsConfig.length()==atomIds.size());
		HashMap<Integer,Integer> chiralMap = new HashMap<Integer, Integer>();
		for(int i=0;i<atomIds.size();i++){
			chiralMap.put(atomIds.get(i),rsConfig.charAt(i)=='R'?OEAtomStereo.RightHanded:OEAtomStereo.LeftHanded);
		}
		for(Integer atomid:chiralMap.keySet()){
			OEAtomBase atm = newMol.GetAtom(new OEHasAtomIdx(atomid));
			OEAtomBaseVector v = new OEAtomBaseVector();
			for(OEAtomBase a:atm.GetAtoms()){
				v.add(a);
			}
			boolean setStereo = atm.SetStereo(v, OEAtomStereo.Tetrahedral, chiralMap.get(atomid));
			if(!setStereo){
				System.err.println("Failed to assign stereo for atom "+atm.GetIdx());
			}else{
			    OEGroupBaseIter andGroups = newMol.GetGroups(new OEHasGroupType(OEGroupType.MDLAndStereo));
			    for(OEGroupBase group:andGroups){
			        if(group.HasAtom(atm)) {
                        group.DeleteAtom(atm);
                    }
                }

//                OEAtomBaseVector orGroup = new OEAtomBaseVector();
//                orGroup.add(atm);
//                OEGroupBase oeGroupBase = newMol.NewGroup(OEGroupType.MDLOrStereo, orGroup);

				OEAtomBaseVector absGroup = new OEAtomBaseVector();
				absGroup.add(atm);
				OEGroupBase oeGroupBase = newMol.NewGroup(OEGroupType.MDLAbsStereo, absGroup);

            }
		}
        newMol.SetTitle(mol.GetTitle()+"_"+rsConfig);
		return newMol;
	}

	public static String[] getAllLists(String[] elements, int lengthOfList)
	{
		//initialize our returned list with the number of elements calculated above
		String[] allLists = new String[(int)Math.pow(elements.length, lengthOfList)];

		//lists of length 1 are just the original elements
		if(lengthOfList == 1) return elements;
		else
		{
			//the recursion--get all lists of length 3, length 2, all the way up to 1
			String[] allSublists = getAllLists(elements, lengthOfList - 1);

			//append the sublists to each element
			int arrayIndex = 0;

			for(int i = 0; i < elements.length; i++)
			{
				for(int j = 0; j < allSublists.length; j++)
				{
					//add the newly appended combination to the list
					allLists[arrayIndex] = elements[i] + allSublists[j];
					arrayIndex++;
				}
			}

			return allLists;
		}
	}


	private OEUnaryAtomPred extractAndFixCommonPiecePredicate(OEGraphMol fixedMol, OEGraphMol targetMol) {
		OEUnaryAtomPred pred = null;
		Random random = new Random();
		random.setSeed(System.currentTimeMillis());
		oechem.OE3DToAtomStereo(fixedMol);
		oechem.OEPerceiveChiral(targetMol);
		HashMap<Integer, Integer> matchList = getFirstMCS(fixedMol, targetMol, OEMCSType.Approximate);
		double maxX = 0.0;
		double maxY = 0.0;
		double maxZ = 0.0;
		double minX = 0.0;
		double minY = 0.0;
		double minZ = 0.0;
		if (matchList != null) {
			for (OEAtomBaseIter iter0 = fixedMol.GetAtoms(); iter0.hasNext(); ) {
				OEAtomBase atom = iter0.next();
				if (matchList.keySet().contains(atom.GetIdx())) {
					Integer idx = matchList.get(atom.GetIdx());
					OEAtomBase targetAtom = targetMol.GetAtom(new OEHasAtomIdx(idx));

					if (targetAtom.IsChiral()) {
						if (targetAtom.HasStereoSpecified(OEAtomStereo.Tetrahedral)) {
							if (atom.IsChiral()) {
								if (atom.HasStereoSpecified(OEAtomStereo.Tetrahedral)) {
									OEAtomBaseVector v1 = new OEAtomBaseVector();
									for (OEAtomBaseIter iter1 = atom.GetAtoms(); iter1.hasNext(); ) {
										v1.add(iter1.next());
									}

									OEAtomBaseVector v2 = new OEAtomBaseVector();
									for (OEAtomBaseIter iter2 = targetAtom.GetAtoms(); iter2.hasNext(); ) {
										v2.add(iter2.next());
									}

									if (atom.GetStereo(v1, OEAtomStereo.Tetrahedral) != targetAtom.GetStereo(v2, OEAtomStereo.Tetrahedral)) {
										System.out.println("No match!");
										return null;
									}
								} else {
									System.out.println("No match!");
									return null;
								}
							} else {
								System.out.println("No match!");
								return null;
							}
						}
					}
				}

			}

			for (OEAtomBaseIter iter = fixedMol.GetAtoms(); iter.hasNext(); ) {
				OEAtomBase atom = iter.next();
				if (matchList.keySet().contains(atom.GetIdx())) {
					Integer idx = matchList.get(atom.GetIdx());

					double[] coords = new double[3];
					fixedMol.GetCoords(atom, coords);
					OEAtomBase targetAtom = targetMol.GetAtom(new OEHasAtomIdx(idx));
					targetMol.SetCoords(targetAtom, coords);

					if (pred == null) {
						pred = new OEHasAtomIdx(idx);
						maxX = coords[0];
						maxY = coords[1];
						maxZ = coords[2];
						minX = coords[0];
						minY = coords[1];
						minZ = coords[2];
					} else {
						pred = new OEOrAtom(pred, new OEHasAtomIdx(idx));
						if (coords[0] > maxX) {
							maxX = coords[0];
						}
						if (coords[1] > maxY) {
							maxY = coords[1];
						}
						if (coords[2] > maxZ) {
							maxZ = coords[2];
						}

						if (coords[0] < minX) {
							minX = coords[0];
						}
						if (coords[1] < minY) {
							minY = coords[1];
						}
						if (coords[2] < minZ) {
							minZ = coords[2];
						}
					}
				}
			}
		} else {
			return null;
		}

		for (OEAtomBaseIter iter = targetMol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			if (!pred.constCall(atom)) {

				if (atom.IsChiral()) {
					if (atom.HasStereoSpecified(OEAtomStereo.Tetrahedral)) {
						OEAtomBaseVector v1 = new OEAtomBaseVector();
						for (OEAtomBaseIter iter1 = atom.GetAtoms(); iter1.hasNext(); ) {
							v1.add(iter1.next());
						}
						System.out.println(atom.GetStereo(v1, OEAtomStereo.Tetrahedral));
					}
				}


				double[] coords = new double[3];
				coords[0] = 10.0 + minX + random.nextDouble() * (maxX - minX);
				coords[1] = 10.0 + minY + random.nextDouble() * (maxY - minY);
				coords[2] = 10.0 + minZ + random.nextDouble() * (maxZ - minZ);
				targetMol.SetCoords(atom, coords);

				oechem.OE3DToAtomStereo(targetMol);

				if (atom.IsChiral()) {
					if (atom.HasStereoSpecified(OEAtomStereo.Tetrahedral)) {
						OEAtomBaseVector v1 = new OEAtomBaseVector();
						for (OEAtomBaseIter iter1 = atom.GetAtoms(); iter1.hasNext(); ) {
							v1.add(iter1.next());
						}
						System.out.println(atom.GetStereo(v1, OEAtomStereo.Tetrahedral));
					}
				}
			}
		}

		return pred;
	}

	private OEGraphMol extractCommonPiece(OEGraphMol fixedMol, OEGraphMol targetMol) {
		OEGraphMol templateMol = new OEGraphMol(fixedMol);
		HashMap<Integer, Integer> matchList = getFirstMCS(fixedMol, targetMol, OEMCSType.Exhaustive);
		if (matchList != null) {
			for (OEAtomBaseIter iter = templateMol.GetAtoms(); iter.hasNext(); ) {
				OEAtomBase atom = iter.next();
				if (!matchList.keySet().contains(atom.GetIdx())) {
					templateMol.DeleteAtom(atom);
				}
			}

			for (OEBondBaseIter iter = templateMol.GetBonds(); iter.hasNext(); ) {
				OEBondBase bond = iter.next();
				int bgnIdx = bond.GetBgnIdx();
				int endIdx = bond.GetEndIdx();
				if (!matchList.keySet().contains(bgnIdx) || !matchList.keySet().contains(endIdx)) {
					templateMol.DeleteBond(bond);
				}
			}

			oechem.OEAssignAromaticFlags(templateMol);

			return templateMol;
		} else {
			return null;
		}

	}
	/*
	public Molecule minimizeMorphingMol(Molecule patternMol, Molecule originalMol, String pdbString, ProgressReporter progressReporter) {
		if (patternMol == null || originalMol == null) {
			return null;
		}
		OEGraphMol originalOE = convertChemAxonMol(originalMol);
		oechem.OESuppressHydrogens(originalOE);
		OEGraphMol targetOE = convertChemAxonMol(patternMol);
		oechem.OESuppressHydrogens(targetOE);
		if (progressReporter != null) {
			progressReporter.reportProgress("Finding common piece ...",20);
		}

		OEUnaryAtomPred pred = null;
		HashMap<Integer, Integer> matchList = getFirstMCS(originalOE, targetOE, OEMCSType.Exhaustive);
		if (matchList != null) {
			for (OEAtomBaseIter iter = originalOE.GetAtoms(); iter.hasNext(); ) {
				OEAtomBase atom = iter.next();
				if (matchList.keySet().contains(atom.GetIdx())) {
					Integer idx = matchList.get(atom.GetIdx());
					if (pred == null) {
						pred = new OEHasAtomIdx(idx);
					} else {
						pred = new OEOrAtom(pred, new OEHasAtomIdx(idx));
					}
				}
			}
		} else {
			return null;
		}

		if (progressReporter != null) {
			progressReporter.reportProgress("Minimizing within receptor ...", 40);
		}
		OEGraphMol resultMol = minimizeOEMol(targetOE, pred, pdbString);
		if (progressReporter != null) {
			progressReporter.reportProgress("finished.", 100);
		}

		if (resultMol != null) {
			return convertOEChemMol(resultMol);
		} else {
			return null;
		}

	}
	 */

	/*
	public Molecule minimizeMoleculeConstrained(Molecule originalMol, Vector<Integer> matchedList, String pdbString, ProgressReporter progressReporter) {
		if (originalMol == null || matchedList == null || matchedList.size() == 0) {
			return null;
		}
		OEGraphMol originalOE = convertChemAxonMol(originalMol);
		if (progressReporter != null) {
			progressReporter.reportProgress("Using selection as template ...",20);
		}

		OEUnaryAtomPred pred = null;
		for (Integer idx : matchedList) {
			if (pred == null) {
				pred = new OEHasAtomIdx(idx - 1);
			} else {
				pred = new OEOrAtom(pred, new OEHasAtomIdx(idx - 1));
			}
		}

		if (progressReporter != null) {
			progressReporter.reportProgress("Minimizing within receptor ...",40);
		}
		OEGraphMol resultMol = minimizeOEMol(originalOE, pred, pdbString);
		if (progressReporter != null) {
			progressReporter.reportProgress("finished.",100);
		}

		if (resultMol != null) {
			return convertOEChemMol(resultMol);
		} else {
			return null;
		}

	}
	 */

	private void makeCliques(OESurface surf, OEGraphMol mol) {
		HashMap<Integer, OEAtomBase> cache = new HashMap();
		for (OEAtomBaseIter iter = mol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			Integer idx = atom.GetIdx();
			cache.put(idx, atom);
		}
		for (int i = 0; i < surf.GetNumVertices(); i++) {
			OEAtomBase atom = cache.get(surf.GetAtomsElement(i));
			if (atom.IsPolar()) {
				surf.SetVertexCliqueElement(i, 1);
			} else {
				surf.SetVertexCliqueElement(i, 2);
			}
		}
		return;
	}

	private void makeCliquesAlternative(OESurface surf, OEGraphMol mol) {
		HashMap<Integer, OEAtomBase> cache = new HashMap();
		for (OEAtomBaseIter iter = mol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			Integer idx = atom.GetIdx();
			cache.put(idx, atom);
		}
		for (int i = 0; i < surf.GetNumVertices(); i++) {
			OEAtomBase atom = cache.get(surf.GetAtomsElement(i));
			if (atom.IsOxygen() || atom.IsNitrogen() || atom.IsPolarHydrogen()) {
				surf.SetVertexCliqueElement(i, 1);
			} else {
				surf.SetVertexCliqueElement(i, 2);
			}
		}
		return;
	}

	public synchronized Molecule calculateAtomicFavorableSurfaceArea(Molecule ligand, String pdbString) {
		if (ligand != null && pdbString != null) {
			Molecule newMol = ligand.cloneMolecule();
			MolHandler mh = new MolHandler(newMol);
			mh.removeHydrogens();
			OEGraphMol oelig = convertChemAxonMol(newMol);
			OEGraphMol oeprot = new OEGraphMol();
			oemolistream pdbIn = new oemolistream();
			pdbIn.SetFormat(OEFormat.PDB);
			pdbIn.openstring(pdbString);
			oechem.OEReadMolecule(pdbIn, oeprot);
			pdbIn.close();
//            oechem.OESuppressHydrogens(oeprot);

			if (oeprot != null && oelig != null) {
				oechem.OEAssignBondiVdWRadii(oeprot);
				oechem.OEAssignBondiVdWRadii(oelig);
				float[] ligandBuriedSurface = buriedSurfaceAreaPortion(oeprot, oelig);
				float[] protBuriedSurface = buriedSurfaceAreaPortion(oelig, oeprot);
				HashMap<Integer, Double> surfaceMap = calculateFavorableContact(oelig, oeprot, protBuriedSurface, ligandBuriedSurface);
				for (MolAtom atom : newMol.getAtomArray()) {
					if (atom.getAtno() > 1) {
						int index = newMol.indexOf(atom);
						if (surfaceMap.containsKey(index)) {
							double surface = surfaceMap.get(index);
							atom.setExtraLabel(new DecimalFormat("##.#").format(Math.abs(surface)));
							if (surface > 0) {
								atom.setExtraLabelColor(Color.GREEN.getRGB());
							} else {
								atom.setExtraLabelColor(Color.RED.getRGB());
							}
						}
					}
				}
				return newMol;
			}
		}
		return null;
	}

	private float[] calculateAtomicSurfaceArea(OESurface surface, OEGraphMol mol) {
		OEFloatArray areas = new OEFloatArray(surface.GetNumTriangles());
		oespicoli.OECalculateTriangleAreas(surface, areas);

		float[] atomic_areas = new float[mol.GetMaxAtomIdx()];
		for (int i = 0; i < atomic_areas.length; i++) {
			atomic_areas[i] = 0.0f;
		}

		for (int i = 0; i < surface.GetNumTriangles(); i++) {
			int v1 = surface.GetTrianglesElement(i * 3);
			int v2 = surface.GetTrianglesElement(i * 3 + 1);
			int v3 = surface.GetTrianglesElement(i * 3 + 2);

			int a1 = surface.GetAtomsElement(v1);
			int a2 = surface.GetAtomsElement(v2);
			int a3 = surface.GetAtomsElement(v3);

			atomic_areas[a1] += areas.getItem(i) / 3.0;
			atomic_areas[a2] += areas.getItem(i) / 3.0;
			atomic_areas[a3] += areas.getItem(i) / 3.0;
		}
		return atomic_areas;
	}

	//
	//todo: calculate complimentary buried surface area
	//
	private HashMap<Integer, Double> calculateFavorableContact(OEGraphMol ligand, OEGraphMol protein, float[] protBuriedSurface, float[] ligandBuriedSurface) {
		OENearestNbrs nn = new OENearestNbrs(ligand, 8.0);
		OENbrsIter nbrsIter = nn.GetNbrs(protein);
		HashMap<Integer, Vector<Integer>> neighborHash = new HashMap<Integer, Vector<Integer>>();
		for (OENbrs neighbor : nbrsIter) {
			int ligAtmIdx = neighbor.GetBgn().GetIdx();
			int proAtmIdx = neighbor.GetEnd().GetIdx();
			double distance = Math.sqrt(neighbor.GetDist2());
			double distance_cutoff = neighbor.GetBgn().GetRadius() + neighbor.GetEnd().GetRadius() + 1.4;
			if (ligandBuriedSurface[ligAtmIdx] > 0.1 && protBuriedSurface[proAtmIdx] > 0.0 && distance < distance_cutoff) {
				if (neighborHash.containsKey(ligAtmIdx)) {
					if (!neighborHash.get(ligAtmIdx).contains(proAtmIdx)) {
						neighborHash.get(ligAtmIdx).add(proAtmIdx);
					}
				} else {
					neighborHash.put(ligAtmIdx, new Vector<Integer>());
					neighborHash.get(ligAtmIdx).add(proAtmIdx);
				}
			}
		}

		HashMap<Integer, Double> favorableSurfaceArea = new HashMap<Integer, Double>();
		for (OEAtomBase ligAtom : ligand.GetAtoms()) {
			int idx = ligAtom.GetIdx();
			double surface = ligandBuriedSurface[idx];
			if (ligandBuriedSurface[idx] > 0.1 && neighborHash.containsKey(idx)) {
				Vector<Integer> nnList = neighborHash.get(idx);
				int favorCount = 0;
				for (int protAtomIdx : nnList) {
					OEAtomBase protAtom = protein.GetAtom(new OEHasAtomIdx(protAtomIdx));
					favorCount += isContactFavorable(ligAtom, protAtom, oechem.OEGetDistance(ligand, ligAtom, protein, protAtom), ligandBuriedSurface[idx]);
				}
				if (favorCount < 0) {
					surface = surface * -1.0f;
					favorableSurfaceArea.put(idx, surface);
				}
				if (favorCount > 0) {
					favorableSurfaceArea.put(idx, surface);
				}
			}
		}
		return favorableSurfaceArea;
	}

	private boolean isPolarAtom(OEAtomBase atom) {
		if (atom.IsSulfur()) {
			return false;
		}
		if (atom.IsCarbon()) {
			return false;
		}
		if (atom.IsHalogen()) {
			return false;
		}
		if (atom.IsHydrogen()) {
			if (!atom.IsPolarHydrogen()) {
				return false;
			}
		}
		return true;
	}

	private boolean isHBondAcceptor(OEAtomBase atom) {
		if (isPolarAtom(atom) && !atom.IsHydrogen()) {
			if (atom.GetExplicitHCount() > 0 || atom.GetImplicitHCount() > 0) {
				if (atom.GetAtomicNum() == OEElemNo.O) {
					return true;
				}
				return false;
			} else {
				return true;
			}
		} else {
			if (atom.GetAtomicNum() == OEElemNo.F) {
				return true;
			}
		}
		return false;
	}

	private boolean isHBondDonor(OEAtomBase atom) {
		if (isPolarAtom(atom)) {
			if (atom.GetExplicitHCount() > 0) {
				return true;
			}
			if (atom.GetImplicitHCount() > 0) {
				return true;
			}
			if (atom.IsPolarHydrogen()) {
				return true;
			}
		}
		return false;
	}

	private int isContactFavorable(OEAtomBase ligAtm, OEAtomBase protAtm, double distance, double ligSurface) {
		//
		if (isHBondAcceptor(ligAtm) && isHBondDonor(protAtm)) {
			if (protAtm.IsHydrogen()) {
				if (distance < 2.6) {
					return 100;
				} else {
					return 1;
				}
			} else {
				if (distance < 3.6) {
					return 100;
				}
			}
		} else if (isHBondAcceptor(protAtm) && isHBondDonor(ligAtm)) {
			if (ligAtm.IsHydrogen()) {
				if (distance < 2.6) {
					return 100;
				} else {
					return 1;
				}
			} else {
				if (distance < 3.6) {
					return 100;
				} else {
					return 1;
				}
			}
		}

		if (ligSurface > 0.3) {
			if (!isPolarAtom(ligAtm) && !isPolarAtom(protAtm)) {
				return 1;
			} else {
				return -1;
			}
		}

		return 0;
	}

	private float[] buriedSurfaceAreaPortion(OEGraphMol refMol, OEGraphMol buriedMol) {
		OEGraphMol new_receptor = new OEGraphMol(refMol);
		OEGraphMol new_ligand = new OEGraphMol(buriedMol);

		float[] atomic_buried_area = new float[buriedMol.GetMaxAtomIdx()];

		OESurface ligand_surface_before = new OESurface();
		oespicoli.OEMakeMolecularSurface(ligand_surface_before, buriedMol, 0.5f, 1.4f);
		float[] atomic_areas_before = calculateAtomicSurfaceArea(ligand_surface_before, buriedMol);

//        OESurface receptor_surface_before = new OESurface();
//        oespicoli.OEMakeMolecularSurface(receptor_surface_before, refMol);
//
		oechem.OEAddMols(new_receptor, new_ligand);
		OESurface complex_surface = new OESurface();
		oespicoli.OEMakeMolecularSurface(complex_surface, new_receptor, 0.5f, 1.4f);
		float[] atomic_areas_after = calculateAtomicSurfaceArea(complex_surface, new_receptor);

		for (int i = 0; i < atomic_areas_before.length; i++) {
			double r = buriedMol.GetAtom(new OEHasAtomIdx(i)).GetRadius();
			atomic_buried_area[i] = (atomic_areas_before[i] - atomic_areas_after[i + refMol.GetMaxAtomIdx()]) / (float) (4.0 * Math.PI * r * r);
		}

		return atomic_buried_area;
	}


	//calculate atomic buried surface area with the context of refMol

	private float[] buriedSurfaceArea(OEGraphMol refMol, OEGraphMol buriedMol) {
		OEGraphMol new_receptor = new OEGraphMol(refMol);
		OEGraphMol new_ligand = new OEGraphMol(buriedMol);

		float[] atomic_buried_area = new float[buriedMol.GetMaxAtomIdx()];

		OESurface ligand_surface_before = new OESurface();
		oespicoli.OEMakeMolecularSurface(ligand_surface_before, buriedMol, 0.5f, 1.4f);
		float[] atomic_areas_before = calculateAtomicSurfaceArea(ligand_surface_before, buriedMol);

//        OESurface receptor_surface_before = new OESurface();
//        oespicoli.OEMakeMolecularSurface(receptor_surface_before, refMol);
//
		oechem.OEAddMols(new_receptor, new_ligand);
		OESurface complex_surface = new OESurface();
		oespicoli.OEMakeMolecularSurface(complex_surface, new_receptor, 0.5f, 1.4f);
		float[] atomic_areas_after = calculateAtomicSurfaceArea(complex_surface, new_receptor);

		for (int i = 0; i < atomic_areas_before.length; i++) {
			atomic_buried_area[i] = (atomic_areas_before[i] - atomic_areas_after[i + refMol.GetMaxAtomIdx()]);
		}

		return atomic_buried_area;
	}

	public static void testSpicoli() {
		oemolistream lgIn = new oemolistream();
		lgIn.SetFormat(OEFormat.SDF);
		lgIn.open("/Users/jfeng1/ligand.sdf");
		OEGraphMol ligand = new OEGraphMol();
		oechem.OEReadMolecule(lgIn, ligand);
		lgIn.close();

		System.out.println(String.format("%d %d", ligand.NumAtoms(), ligand.GetMaxAtomIdx()));

		oemolistream pdbIn = new oemolistream();
		pdbIn.SetFormat(OEFormat.PDB);
		pdbIn.open("/Users/jfeng1/receptor.pdb");
		OEGraphMol receptor = new OEGraphMol();
		oechem.OEReadMolecule(pdbIn, receptor);
		oechem.OEAssignHybridization(receptor);
		pdbIn.close();

		OEChemFunc.getInstance().buriedSurfaceArea(ligand, receptor);
	}

	public void calculate3DSurfaceVolumeDipole(PropertyMolecule originalMol) throws Exception {
		openeye.oechem.OEStopwatch stopwatch = new OEStopwatch();
		if (originalMol != null) {
			System.out.println("Starting volume calculation.");
			stopwatch.Start();
			OEGraphMol oemol = originalMol.getMol3d();
			if(oemol!=null) {
				oechem.OEAddExplicitHydrogens(oemol);
				oechem.OEAssignBondiVdWRadii(oemol);
				OESurface surface = new OESurface();
				oespicoli.OEMakeAccessibleSurface(surface, oemol, 0.5f, 1.4f);
				//oespicoli.OEMakeMolecularSurface(surface, oemol, 0.5f, 1.4f);
				//makeCliquesAlternative(surface,oemol);
				makeCliquesAlternative(surface, oemol);

				float volume = oespicoli.OESurfaceVolume(surface);
				float totalSurfaceArea = oespicoli.OESurfaceArea(surface);
				float polarSurfaceArea = oespicoli.OESurfaceCliqueArea(surface, 1);
//				float hydrophobicSurfaceArea = oespicoli.OESurfaceCliqueArea(surface, 2);
				double dipole_moment = calculateDipoleMoment(oemol);

				NumberFormat nf = new DecimalFormat("##.###");

				//"Volume3D", "SurfaceArea3D","PolarSurfaceArea3D", "DipoleMoment3D"
				originalMol.addProperty("Volume3D", nf.format(volume));
				originalMol.addProperty("SurfaceArea3D", nf.format(totalSurfaceArea));
				originalMol.addProperty("PolarSurfaceArea3D", nf.format(polarSurfaceArea));
//				originalMol.addProperty("PolarSurfacePercentage", nf.format(polarSurfaceArea * 100 / totalSurfaceArea));
//				originalMol.addProperty("HydrophobicSurfaceArea3D", nf.format(hydrophobicSurfaceArea));
				originalMol.addProperty("DipoleMoment3D", nf.format(dipole_moment));
				System.out.println(stopwatch.Elapsed() + " elapsed");
			}
		}
	}

	public Molecule calculateQuacpacChargeMMFF94(Molecule mol) {
		return calculateQuacpacCharge(mol, OECharges.MMFF94);
	}

	public Molecule calculateQuacpacChargeAM1BCC(Molecule mol) {
		return calculateQuacpacCharge(mol, OECharges.AM1BCCSPt);
	}

	public double calculateDipoleMoment(OEGraphMol oemol) {
		double dipoleMoment = 0.0;
		if (oemol != null) {
			double ux = 0.0;
			double uy = 0.0;
			double uz = 0.0;
			boolean noHydrogen = false;
			boolean debug = false;
			oequacpac.OEAssignPartialCharges(oemol, OECharges.AM1BCCSPt, noHydrogen, debug);
			double[] baseXyz = new double[3 * oemol.NumAtoms()];
			oemol.GetCoords(baseXyz);
			double cx = 0.0;
			double cy = 0.0;
			double cz = 0.0;
			for (int i = 0; i < oemol.NumAtoms(); i++) {
				cx += baseXyz[i * 3];
				cy += baseXyz[i * 3 + 1];
				cz += baseXyz[i * 3 + 2];
			}
			cx /= oemol.NumAtoms();
			cy /= oemol.NumAtoms();
			cz /= oemol.NumAtoms();
			OEAtomBaseIter iterator = oemol.GetAtoms();
			while (iterator.hasNext()) {
				OEAtomBase oeAtom = iterator.next();
				double charge = oeAtom.GetPartialCharge();
				double[] xyz = new double[3];
				oemol.GetCoords(oeAtom, xyz);
//                System.out.println(String.format("%s %f %f %f %f",oeAtom.GetAtomicNum(),xyz[0]-cx,xyz[1]-cy,xyz[2]-cz,charge));
				ux += charge * (xyz[0] - cx);
				uy += charge * (xyz[1] - cy);
				uz += charge * (xyz[2] - cz);
			}
			dipoleMoment = Math.sqrt(ux * ux + uy * uy + uz * uz) / (0.20822678);   //convert to Debye
		}
		return dipoleMoment;
	}

	public synchronized Molecule calculateQuacpacCharge(Molecule mol, int charge_type) {
		if (mol != null) {
			Molecule newMol = mol.cloneMolecule();
			OEGraphMol oemol = convertChemAxonMol(newMol);
			if (oemol != null) {
				HashMap<Integer, Double> chargeMap = getOEChargeMap(charge_type, oemol);
				for (MolAtom atom : newMol.getAtomArray()) {
					if (atom.getAtno() > 1) {
						int index = newMol.indexOf(atom);
						if (chargeMap.containsKey(index)) {
							double charge = chargeMap.get(index);
							atom.setExtraLabel(new DecimalFormat("##.##").format(charge));
							if (charge > 0) {
								atom.setExtraLabelColor(Color.BLUE.getRGB());
							} else {
								atom.setExtraLabelColor(Color.RED.getRGB());
							}
						}
					}
				}
				return newMol;
			}
		}
		return null;
	}

//	public String generate_iupac_name(Molecule mol) {
//		OEGraphMol oemol = convertChemAxonMol(mol);
//		String name = oeiupac.OECreateIUPACName(oemol);
//		if (name.contains("BLAH")) {
//			return "ERROR";
//		} else {
//			return name;
//		}
//	}

	public HashMap<Integer, Double> getOEChargeMap(int charge_type, OEGraphMol oemol) {
		HashMap<Integer, Double> chargeMap = new HashMap<Integer, Double>();
		boolean noHydrogen = true;
		boolean debug = false;
		oequacpac.OEAssignPartialCharges(oemol, charge_type, noHydrogen, debug);
		OEAtomBaseIter iterator = oemol.GetAtoms();
		while (iterator.hasNext()) {
			OEAtomBase oeAtom = iterator.next();
			System.out.println(String.format("%d %f", oeAtom.GetIdx(), oeAtom.GetPartialCharge()));
			chargeMap.put(oeAtom.GetIdx(), oeAtom.GetPartialCharge());
		}
		return chargeMap;
	}

	/*
"""
Molecular Surface Properties.

Compute the following properties for a input oemol()

Surface Area
Polar Surface Area
Volume
Hydrophobic Area


"""
from openeye.oechem import *
from openeye.oespicoli import *
import sys
def makeCliques( surf, mol):
    """
    Hydrophillic is 1
    Hydrophobic is 2
    Using the oechem definitions of polar, which is
    any atom that is not a C or a H
    """
    cache = dict([(atom.GetIdx(), atom) for atom in mol.GetAtoms()])
    for i in xrange(surf.GetNumVertices()):
        atom = cache[surf.GetAtomsElement(i)]
        if atom.IsPolar():
            surf.SetVertexCliqueElement(i, 1)
#                surf.SetColorElement(i, 0.0, 0.0, 1.0)
        else:
#                surf.SetColorElement(i, 1.0, 0.0, 0.0)
            surf.SetVertexCliqueElement(i,2)
    return True


def makeCliquesAlternative( surf, mol):
    """
    Hydrophillic is 1
    Hydrophobic is 2
    Using the oechem definitions of polar, which is
    any atom that is not a C or a H
    """
    cache = dict([(atom.GetIdx(), atom) for atom in mol.GetAtoms()])
    for i in xrange(surf.GetNumVertices()):
        atom = cache[surf.GetAtomsElement(i)]
        if atom.IsOxygen() or atom.IsNitrogen() or atom.IsPolarHydrogen():
            surf.SetVertexCliqueElement(i, 1)
#                surf.SetColorElement(i, 0.0, 0.0, 1.0)
        else:
#                surf.SetColorElement(i, 1.0, 0.0, 0.0)
            surf.SetVertexCliqueElement(i,2)
    return True



def main(argv=[__name__]):
    """

    """
    OEAssignBondiVdWRadii(mol)
    surface = OESurface()
    #make a surface with a 1.4A probe
    OEMakeMolecularSurface(surface, mol, 0.5, 1.4)
    #we have 2 clique methods
    #makeCliques(surface,mol)
    makeCliquesAlternative(surface,mol)

    volume = OESurfaceVolume(surface)
    totalSurfaceArea = OESurfaceArea(surface)
    polarSurfaceArea = OESurfaceCliqueArea(surface,1)
    hydrophibicSurfaceArea = OESurfaceCliqueArea(surface, 2)
    #could also simply do hydrophbic = total - polar

    print('Volume %8.3f TotalSurfaceArea %8.3f PolarSurfaceArea %8.3f HydrophobicSurfaceArea %8.3f' % (volume,\
                                                                                                       totalSurfaceArea,\
                                                                                                       polarSurfaceArea,\
                                                                                                       hydrophibicSurfaceArea))
    return






if __name__ == '__main__':
    mol = OEMol()
    ifs = oemolistream(sys.argv[1])
    OEReadMolecule(ifs, mol)
    ifs.close()
    main(mol)


     */

	/*
    public Molecule getMorphingMol2(Molecule patternMol, Molecule originalMol, ProgressReporter progressReporter) {
        if (patternMol == null || originalMol == null) {
            return null;
        }

        OEGraphMol originalOE = convertChemAxonMol(originalMol);
        oechem.OE3DToAtomStereo(originalOE);
        OEGraphMol targetOE = convertChemAxonMol(patternMol);
        oechem.OEPerceiveChiral(targetOE);
        if(progressReporter!=null){
            progressReporter.reportProgress(10,"Finding common piece ...");
        }
        OEGraphMol templateOE = extractCommonPiece(originalOE,targetOE);
        Molecule bestMol = null;
        OEGraphMol oeBestMol = null;
        double bestScore = 0.0;
        if(templateOE!=null){
            OEMol oeMol = new OEMol(targetOE);
            if(progressReporter!=null){
                progressReporter.reportProgress(10,"Generating restrictive conformers ...");
            }
            getRestrictiveConformers(oeMol,templateOE,50,"mmffs");
            if(oeMol.NumAtoms()>0){
                for (OEConfBaseIter iter = oeMol.GetConfs(); iter.hasNext();) {
                    OEGraphMol oemol = new OEGraphMol(iter.next());
                    double score = getComboScoreOE(oemol,originalOE);
                    if(score > bestScore){
                        oeBestMol = oemol;
                        bestScore = score;
                    }
                }
                if(oeBestMol!=null){
                    bestMol = convertOEChemMol(oeBestMol);
                }
            }else{
                if(progressReporter!=null){
                    progressReporter.reportProgress(20,"Failed to find matching piece, running shape overlay ...");
                }
                Vector<SuperpositionSolution> solutions = colorOverlapRigid2Flex(templateOE, targetOE, 50, progressReporter, true);
                if(solutions!=null&&solutions.size()>0){
                    SuperpositionSolution solution = solutions.get(solutions.size()-1);
                    bestMol = solution.getTarget();
                }
            }
        }



//        Vector<SuperpositionSolution> solutions = colorOverlapRigid2Flex(templateOE, targetOE, 50, progressReporter, true);
//        if(solutions!=null&&solutions.size()>0){
//            SuperpositionSolution solution = solutions.get(solutions.size()-1);
//            targetMol = solution.getTarget();
//        }

        return bestMol;

    }
*/



	public static void main(String[] args) {
		OEChemFunc.initialize();
		String smiles = "CCC(=O)N[C@@H](CCC(=O)[O-])C(=O)N[C@@]1(CCCCC[C@](NC(=O)[C@@H](NC1=O)CCCNC(=[NH2+])N)(Cc2c[nH]cn2)C(=O)N[C@@H](CC(=O)[O-])C=O)Cc3ccccc3";
		String smiles2 = "CCC(=O)N[C@@H](CCC([O-])=O)C(=O)N[C@]1(CC2=CC=CC=C2)CCCCCCNC(=O)[C@H](CC([O-])=O)NC(=O)[C@H](CC2=CNC=N2)NC(=O)[C@H](CCCNC(N)=[NH2+])NC1=O";
		OEGraphMol mol = new OEGraphMol();
		oechem.OEParseSmiles(mol,smiles2);
        PropertyMolecule m = new PropertyMolecule(mol);
		Vector<OEGraphMol> molList = getChiralMols(m.getMol(),false);
		for(OEGraphMol m1:molList){
			String smi = oechem.OEMolToSmiles(m1);
			System.out.println(smi);
		}


//        testSpicoli();
//        try {
//            Molecule targetMol = MolImporter.importMol(FileUtils.readFileToString(new File("/Users/feng/test.sdf")));
//            Molecule mol = OEChemFunc.getInstance().calculateQuacpacCharge(targetMol, OECharges.AM1BCCSPt);
//        } catch (IOException e) {
//            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }

//        testMCS();
		//testPred();
		//testFit();
		/*
        Molecule targetMol = null;
        Molecule templateMol = null;
        try {                                                                                                                              g
//            targetMol = MolImporter.importMol(FileUtils.readFileToString(new File("/Users/feng/tmp/846612_2D.sdf")));
//            templateMol = MolImporter.importMol(FileUtils.readFileToString(new File("/Users/feng/tmp/846612_xray.sdf")));
            targetMol = MolImporter.importMol(FileUtils.readFileToString(new File("/Users/feng/tmp/840586_hack3D.sdf")));
            templateMol = MolImporter.importMol(FileUtils.readFileToString(new File("/Users/feng/tmp/840586_template.sdf")));
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
        OEGraphMol template = OEChemFunc.getInstance().convertChemAxonMol(templateMol);
        OEGraphMol target = OEChemFunc.getInstance().convertChemAxonMol(targetMol);
        OEMol refMol = new OEMol(template);
        OEMol fitMol = new OEMol(target);
        */
		/*
        try {
            System.out.println("reading file ...");
            String pdbString = FileUtils.readFileToString(new File("/Users/feng/PB2_receptor.pdb"));
            byte[] result = OEChemFunc.getInstance().calcElectrostaticMap(pdbString, OEFormat.PDB, true);
//            FileUtils.writeByteArrayToFile(new File("junk.phi"), result );
        } catch (IOException e) {
            e.printStackTrace();
        }
        */
		//testOMega();
		/*
        Vector<SuperpositionSolution> results = oechemFunc.colorOverlapRigid2Flex(template,target,3000,new ProgressReporter() {
            public void reportProgress(int i, Object message) {
                System.out.println(message);
            }

            public boolean isCanceled() {
                return false;
            }
        });

        for (int i = 0; i < results.size(); i++) {
            SuperpositionSolution solution = results.elementAt(i);
            System.out.println(solution.getShapeScore());
        }*/
	}

	private static void testPred() {
		oemolistream ifs = new oemolistream();
		ifs.open("/Users/feng/test.sdf");
		OEMol mol = new OEMol();
		oechem.OEReadMolecule(ifs, mol);
		OEGraphMol oegraphMol = new OEGraphMol(mol);
		Vector<Integer> idList = new Vector<>();
		idList.add(0);
		idList.add(8);
//		OEUnaryAtomPred pred = new OEHasAtomIdx(0);
//		pred = new OEOrAtom(pred, new OEHasAtomIdx(8));
		OEGraphMol minimizedMol = OEChemFunc.getInstance().minimizeOEMol(oegraphMol, idList);
		System.out.println(OEChemFunc.getInstance().convertOEChemMol(minimizedMol).toFormat("sdf"));

	}

	private static void testLowestEnergy() {
		oemolistream ifs = new oemolistream();
		ifs.open("/Users/feng/test.sdf");
		OEMol mol = new OEMol();
		int id = 0;
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileWriter("/Users/feng/confs.out"));
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}

		OEChemFunc oechemFunc = OEChemFunc.getInstance();
		while (oechem.OEReadMolecule(ifs, mol)) {
			double low_2000 = oechemFunc.generateMultipleConformers(mol, 2000,100);
			double low_40 = oechemFunc.generateMultipleConformers(mol, 40,100);
			out.println(String.format("%d %f %f", id++, low_40, low_2000));
		}

		out.close();
	}

	private static void testMCS() {
		oemolistream ifs = new oemolistream();
		ifs.open("/Users/jfeng/template.mol");
		OEGraphMol template = new OEGraphMol();
		oechem.OEReadMolecule(ifs, template);
		System.out.println(template.NumAtoms());
		ifs.close();

		ifs.open("/Users/jfeng/target.mol");
		OEGraphMol target = new OEGraphMol();
		oechem.OEReadMolecule(ifs, target);
		oechem.OESuppressHydrogens(target);
		System.out.println(target.NumAtoms());
		ifs.close();

		OEChemFunc.getInstance().getFirstMCS(template, target, OEMCSType.Exhaustive);

	}

	private static void testFit() {
		oemolistream ifs = new oemolistream();
		oemolistream ifs2 = new oemolistream();
//        ifs2.open("/Users/feng/tmp/846612_xray.sdf");
		ifs2.open("/Users/jfeng1/template.sdf");
		OEMol refMol = new OEMol();
		oechem.OEReadMolecule(ifs2, refMol);

//        ifs.open("/Users/feng/tmp/846612_2D.sdf");
		ifs.open("/Users/jfeng1/target2.smi");

		OEMol fitMol = new OEMol();
		oechem.OEReadMolecule(ifs, fitMol);

		OEGraphMol piece = OEChemFunc.getInstance().extractCommonPiece(new OEGraphMol(refMol), new OEGraphMol(fitMol));

		OEOmega omega = new OEOmega();

		omega.GetOptions().SetFixMol(piece);

		try {
			omega.call(fitMol);
			Molecule mol = OEChemFunc.getInstance().convertOEChemMol(new OEGraphMol(fitMol));
			System.out.println(fitMol.NumAtoms() + " found.");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	private static void testOMega() {
		oemolistream ifs = new oemolistream();
		oemolistream ifs2 = new oemolistream();
		ifs2.open("/Users/jfeng1/tmp/template.sdf");
		OEMol refMol = new OEMol();
		oechem.OEReadMolecule(ifs2, refMol);

		ifs.open("/Users/jfeng1/tmp/hack2D.sdf");
		OEMol fitMol = new OEMol();
		oechem.OEReadMolecule(ifs, fitMol);


		OEOmega omega = new OEOmega();
		omega.GetOptions().SetFromCT(true);
		omega.GetOptions().SetMaxConfs(2000);
		omega.GetOptions().SetEnergyWindow(10);
		omega.GetOptions().SetBuildForceField(OEForceFieldType.MMFF94S);
		omega.GetOptions().SetSearchForceField(OEForceFieldType.MMFF94S);

		OEBestOverlay best = new OEBestOverlay();
		best.SetRefMol(refMol);
		best.SetColorForceField(OEColorFFType.ImplicitMillsDean);
		best.SetColorOptimize(true);
		System.out.println("Self color:" + best.GetRefSelfColor());

		omega.call(fitMol);

		OEBestOverlayResultsIter resiter = best.Overlay(fitMol);
		OEBestOverlayScoreIter scoreiter = new OEBestOverlayScoreIter();
		oeshape.OESortOverlayScores(scoreiter, resiter, new OEHighestComboScore());
		OEBestOverlayScore score;
		OEGraphMol outMol = null;
		int count = 0;
		for (; scoreiter.IsValid(); scoreiter.Increment()) {
			score = scoreiter.Target();
			outMol = new OEGraphMol(fitMol.GetConf(new OEHasConfIdx(score.getFitconfidx())));
			score.Transform(outMol);
			System.out.println(String.format("%d %.4f", count, score.GetComboScore() / 2.0));
			if (count++ == 5) {
				break;
			}
		}
	}

}

