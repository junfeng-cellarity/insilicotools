package com.insilico.application.insilicotools.util;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.data.pKa;
import com.insilico.application.insilicotools.data.pKa_type;
import openeye.oechem.*;

import javax.vecmath.Vector3d;
import java.util.HashMap;
import java.util.Vector;

public class Protonator {
	OESubSearch aniline_pattern;
	Vector<ChargePattern> charget_list;
	static Protonator _this;

	public Protonator() {
		aniline_pattern = new OESubSearch();
		aniline_pattern.Init("[NH2]-*(=,:*)-,:*");
//        String chargedAtoms = "[#7;!$([#7]=,:,#*);!$([#7]-,:*=,:*);!$(N-N)] +1 basic amines\n" +
//                "[N;$([ND1]=C-N)] +1 amidenes\n" +
//                "[#8;$([OD1]-[C,S,P]=O)] -1 acids";
		charget_list = new Vector<ChargePattern>() {
			{
				add(new ChargePattern("[#7;!$([#7]=,:,#*);!$([#7]-,:*=,:*);!$(N-N)]", 1));
				add(new ChargePattern("[N;$([ND1]=C-N)]", 1));
				add(new ChargePattern("[#8;$([OD1]-[C,S,P]=O)]", -1));
			}
		};
	}

	public static Protonator getInstance() {
		if (_this == null) {
			_this = new Protonator();
		}
		return _this;
	}

	public boolean adjust_hydrogens(OEGraphMol mol, boolean protonate) {
		boolean adjusted = false;
		if (protonate) {
			Vector<PropertyMolecule> molList = new Vector<PropertyMolecule>();
			PropertyMolecule pmol = new PropertyMolecule(mol);
			molList.add(pmol);
			ChemFunc.calculateChemAxon_pKa(molList,null);
			HashMap<Integer,pKa> pKaDict = pmol.getpKaMap(PropertyMolecule.CHEMAXON_PKA);

			for (OEAtomBase atom : mol.GetAtoms()) {
				atom.SetFormalCharge(0);
			}
			for (ChargePattern pattern : charget_list) {
				for (OEMatchBase p : pattern.pat.Match(mol)) {
					for (OEMatchPairAtom pair : p.GetAtoms()) {
						OEAtomBase atom = pair.getTarget();
						if(pKaDict.containsKey(atom.GetIdx())){
							pKa pKa = pKaDict.get(atom.GetIdx());
							if(pKa.getType()== pKa_type.BASIC&&pKa.getValue()<7.0){
								continue;
							}
							int h_count = oechem.OEDefaultMDLHCount(atom) - atom.GetExplicitValence() + pattern.charge;
							h_count = h_count < 0 ? 0 : h_count;
							atom.SetImplicitHCount(h_count);
							atom.SetFormalCharge(pattern.charge);
							adjusted = true;
						}
					}
				}
			}
		}
		oechem.OEAssignMDLHydrogens(mol);
		if(mol.GetDimension()==3) {
			oechem.OEAddExplicitHydrogens(mol);
			oechem.OESet3DHydrogenGeom(mol);
		}
		return adjusted;
	}

	public void fix_aniline_hydrogens(OEGraphMol mol) {
		HashMap<Integer, Integer> used = new HashMap<Integer, Integer>();
		for (OEMatchBase p : aniline_pattern.Match(mol)) {
			Vector<OEAtomBase> atomList = new Vector<OEAtomBase>();
			for (OEMatchPairAtom pair : p.GetAtoms()) {
				atomList.add(pair.getTarget());
			}
			OEAtomBase N_Atm = atomList.get(0);
			OEAtomBase C1_Atm = atomList.get(1);
			OEAtomBase C2_Atm = atomList.get(2);
			OEAtomBase C3_Atm = atomList.get(3);
			int N_idx = N_Atm.GetIdx();
			if (!used.containsKey(N_idx)) {
				used.put(N_idx, 1);
				Vector<Vector3d> v = new Vector<Vector3d>();
				double[] n_crd = new double[3];
				mol.GetCoords(N_Atm, n_crd);
				double[] c1_crd = new double[3];
				mol.GetCoords(C1_Atm, c1_crd);
				double[] c2_crd = new double[3];
				mol.GetCoords(C2_Atm, c2_crd);
				double[] c3_crd = new double[3];
				mol.GetCoords(C3_Atm, c3_crd);

				Vector3d N_crd = new Vector3d(n_crd);
				Vector3d C1_crd_a = new Vector3d(c1_crd);
				Vector3d C1_crd_b = new Vector3d(c1_crd);
				Vector3d C2_crd = new Vector3d(c2_crd);
				Vector3d C3_crd = new Vector3d(c3_crd);

				C1_crd_a.sub(C2_crd);
				C1_crd_a.normalize();
				C1_crd_a.scale(1.020);
				v.add(C1_crd_a);

				C1_crd_b.sub(C3_crd);
				C1_crd_b.normalize();
				C1_crd_b.scale(1.020);
				v.add(C1_crd_b);

				int idx = 0;
				for (OEAtomBase nbr : N_Atm.GetAtoms()) {
					if (nbr.GetAtomicNum() == 1) {
						Vector3d H_crd = new Vector3d(N_crd);
						H_crd.add(v.get(idx));
						idx += 1;
					}
				}
			}
		}
	}

	public boolean protonate(OEGraphMol mol) {
		oechem.OESuppressHydrogens(mol);
		boolean adjusted = adjust_hydrogens(mol, true);
		fix_aniline_hydrogens(mol);
		return adjusted;
	}


	class ChargePattern {
		String smarts;
		int charge;
		OESubSearch pat;

		ChargePattern(String smarts, int charge) {
			this.smarts = smarts;
			this.charge = charge;
			this.pat = new OESubSearch();
			pat.Init(smarts);
		}

	}
}
