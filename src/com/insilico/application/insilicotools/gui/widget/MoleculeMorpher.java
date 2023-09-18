package com.insilico.application.insilicotools.gui.widget;

import chemaxon.formats.MolImporter;
import chemaxon.struc.CTransform3D;
import chemaxon.struc.MolAtom;
import chemaxon.struc.MolBond;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.util.Protonator;
import openeye.oechem.*;

import java.util.Vector;

public class MoleculeMorpher {

	public static String linkMolecule(String parentString, String fragmentString, int atomId1, int atomId2, int length) throws Exception {
		return linkMolecule(parentString, fragmentString, atomId1, atomId2, length, false);
	}

	public static String linkMolecule(String parentString, String fragmentString, int atomId1, int atomId2, int length, boolean force) throws Exception {
		OEChemFunc oechemFunc = OEChemFunc.getInstance();
		Molecule parent = MolImporter.importMol(parentString);
		Molecule fragment = MolImporter.importMol(fragmentString);
		if (parent != null && fragment != null) {
			MolAtom headAtom = parent.getAtom(atomId1 - 1);
			MolAtom tailAtom = fragment.getAtom(atomId2 - 1);
			if (headAtom == null || tailAtom == null) {
				throw new Exception("No Atom found.");
			}


			MolAtom attachPoint1 = parent.getAtom(atomId1 - 1);
			if (attachPoint1.getAtno() == 1) {
				throw new Exception("Hydrogen can't be used as attaching point.");
			}

			int explictHCount = 0;
			for (int i = 0; i < attachPoint1.getBondCount(); i++) {
				MolBond bond = attachPoint1.getBond(i);
				MolAtom neighborAtom = bond.getOtherAtom(attachPoint1);
				if (neighborAtom.getAtno() == 1) {
					explictHCount += 1;
				}
			}

			if ((attachPoint1.getImplicitHcount() == 0 && explictHCount == 0) || attachPoint1.getBondCount() == 0) {
				if (!force) {
					throw new Exception("No free valence to add new group.");
				}
			}

			MolAtom attachPoint2 = fragment.getAtom(atomId2 - 1);
			if (attachPoint2.getAtno() == 1) {
				throw new Exception("Hydrogen can't be used as attaching point.");
			}

			explictHCount = 0;
			for (int i = 0; i < attachPoint2.getBondCount(); i++) {
				MolBond bond = attachPoint2.getBond(i);
				MolAtom neighborAtom = bond.getOtherAtom(attachPoint2);
				if (neighborAtom.getAtno() == 1) {
					explictHCount += 1;
				}
			}

			if ((attachPoint2.getImplicitHcount() == 0 && explictHCount == 0) || attachPoint2.getBondCount() == 0) {
				if (!force) {
					throw new Exception("No free valence to add new group.");
				}
			}

			//Remove a hydrogen atom from the parent molecule
			MolAtom removedHatom1 = null;
			for (int i = 0; i < attachPoint1.getBondCount(); i++) {
				MolBond bond = attachPoint1.getBond(i);
				MolAtom neighborAtom = bond.getOtherAtom(attachPoint1);
				if (neighborAtom.getAtno() == 1) {
					parent.removeBond(bond);
					parent.removeAtom(neighborAtom);
					removedHatom1 = neighborAtom;
					break;
				}
			}

			MolAtom removedHatom2 = null;
			for (int i = 0; i < attachPoint2.getBondCount(); i++) {
				MolBond bond = attachPoint2.getBond(i);
				MolAtom neighborAtom = bond.getOtherAtom(attachPoint2);
				if (neighborAtom.getAtno() == 1) {
					fragment.removeBond(bond);
					fragment.removeAtom(neighborAtom);
					removedHatom2 = neighborAtom;
					break;
				}
			}

			//Combine parent molecule with fragment molecule
			for (MolAtom atom : fragment.getAtomArray()) {
				parent.insertAtom(parent.getAtomCount(), atom);
			}

			for (MolBond bond : fragment.getBondArray()) {
				parent.insertBond(parent.getBondCount(), bond);
			}
			if (length == -1) {
				length = getPreferredLength(attachPoint1, attachPoint2);
			}
			if (length > 0) {
				OEUnaryAtomPred pred = null;
				Vector<Integer> idList = new Vector<Integer>();
				MolAtom attachPoint3 = null;
				MolAtom attachPoint4 = null;
				StringBuilder sb = new StringBuilder();
				for(int i=0;i<length+2;i++){
					sb.append("C");
				}
				String fragSmi = sb.toString();
				Molecule fragMol = MolImporter.importMol(fragSmi);
				fragMol.clean(3, null);
				for (MolAtom atom : fragMol.getAtomArray()) {
					if (fragMol.indexOf(atom) == 0) {
						attachPoint3 = atom;
					}
					if (fragMol.indexOf(atom) == length + 1) {
						attachPoint4 = atom;
					}
					parent.insertAtom(parent.getAtomCount(), atom);
					System.out.println(parent.indexOf(atom) + " should be optimized.");
				}

				if (attachPoint3 == null || attachPoint4 == null) {
					throw new Exception("No Attach point defined in the group.");
				}

				if (removedHatom1 != null && removedHatom2 != null) {
					CTransform3D transform3D = new CTransform3D();
					transform3D.setTranslation((removedHatom1.getX() + removedHatom2.getX()) / 2.0, (removedHatom1.getY() + removedHatom2.getY()) / 2.0, (removedHatom1.getZ() + removedHatom2.getZ()) / 2.0);
					fragMol.transform(transform3D);
				}


				if (removedHatom1 != null) {
					attachPoint3.setXYZ(attachPoint1.getX(), attachPoint1.getY(), attachPoint1.getZ());
				}

				if (removedHatom2 != null) {
					attachPoint4.setXYZ(attachPoint2.getX(), attachPoint2.getY(), attachPoint2.getZ());
				}

				int indexOfPoint3 = fragMol.indexOf(attachPoint3);
				int indexOfPoint4 = fragMol.indexOf(attachPoint4);
				OEGraphMol oemol = oechemFunc.convertChemAxonMol(fragMol);
				Vector<Integer> ids = new Vector<>();
				ids.add(indexOfPoint3);
				ids.add(indexOfPoint4);
				oechem.OEAddExplicitHydrogens(oemol);
				Protonator.getInstance().protonate(oemol);
				OEGraphMol minimizedMol = oechemFunc.minimizeOEMol(oemol, ids);
				oechem.OESuppressHydrogens(oemol);
				fragMol = oechemFunc.convertOEChemMol(minimizedMol);
				attachPoint3 = fragMol.getAtom(indexOfPoint3);
				attachPoint4 = fragMol.getAtom(indexOfPoint4);

				MolAtom newPoint3 = null;
				MolAtom newPoint4 = null;
				for (int j = 0; j < attachPoint3.getBondCount(); j++) {
					MolAtom neighbor = attachPoint3.getBond(j).getOtherAtom(attachPoint3);
					if (neighbor.getAtno() > 1) {
						newPoint3 = neighbor;
						break;
					}
				}

				for (int j = 0; j < attachPoint4.getBondCount(); j++) {
					MolAtom neighbor = attachPoint4.getBond(j).getOtherAtom(attachPoint3);
					if (neighbor.getAtno() > 1) {
						newPoint4 = neighbor;
						break;
					}
				}

				if (newPoint3 == null || newPoint4 == null) {
					throw new Exception("Unable to link molecules.");
				}


				/*
								if (removedHatom2 != null) {
									CTransform3D transform3D = new CTransform3D();
									transform3D.setTranslation(removedHatom2.getX(), removedHatom2.getY(), removedHatom2.getZ());
									attachPoint4.transform(transform3D, false);
								}
								*/
				for (MolAtom atom : fragMol.getAtomArray()) {
					if (atom == attachPoint3 || atom == attachPoint4) {
						continue;
					}
					parent.insertAtom(parent.getAtomCount(), atom);
				}

				for (MolBond bond : fragMol.getBondArray()) {
					if (bond.getAtom1() == attachPoint3 || bond.getAtom2() == attachPoint3) {
						continue;
					}
					if (bond.getAtom1() == attachPoint4 || bond.getAtom2() == attachPoint4) {
						continue;
					}
					parent.insertBond(parent.getBondCount(), bond);
				}


				MolBond bond1 = new MolBond(attachPoint1, newPoint3);
				parent.insertBond(parent.getBondCount(), bond1);

				MolBond bond2 = new MolBond(attachPoint2, newPoint4);
				parent.insertBond(parent.getBondCount(), bond2);

				oemol = oechemFunc.convertChemAxonMol(parent);
				oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol);
				oechem.OESuppressHydrogens(oemol);
				oechem.OEAddExplicitHydrogens(oemol);
				Molecule molecule = oechemFunc.convertOEChemMol(oemol);
				return molecule.toFormat("mol");
			} else {
				MolBond bond = new MolBond(attachPoint1, attachPoint2);
				parent.insertBond(parent.getBondCount(), bond);
				OEGraphMol oemol = oechemFunc.convertChemAxonMol(parent);
				oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol);
				oechem.OESuppressHydrogens(oemol);
				oechem.OEAddExplicitHydrogens(oemol);
				Molecule molecule = oechemFunc.convertOEChemMol(oemol);
				return molecule.toFormat("mol");
			}
		}
		return null;
	}

	private static int getPreferredLength(MolAtom attachPoint1, MolAtom attachPoint2) {
		int length;

		if (attachPoint1.isBoundTo(attachPoint2)) {
			return 4;
		}

		double distance = getDistance(attachPoint1, attachPoint2);
		if (distance < 2.0) {
			return 0;
		}
		length = (int) (distance / 1.25) - 1;
		return length;
	}

	public static double getDistance(MolAtom atom1, MolAtom atom2) {
		if (atom1 != null && atom2 != null) {
			double x1 = atom1.getX();
			double y1 = atom1.getY();
			double z1 = atom1.getZ();
			double x2 = atom2.getX();
			double y2 = atom2.getY();
			double z2 = atom2.getZ();
			return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
		}
		return 0.0;
	}

	public static String attachGroupToMolecule(String molString, int atomId, MorphingFragment fragment) throws Exception {
		Molecule mol = MolImporter.importMol(molString);
		Molecule frag;
		if (fragment.getMolFile() == null) {
			frag = MolImporter.importMol(fragment.getSmiles());
			frag.clean(3, null);
//            frag = ModelingTasks.getCorina3DStructure(fragment.getSmiles()).get(0);
			fragment.setMolFile(frag.toFormat("mol"));
		} else {
			frag = MolImporter.importMol(fragment.getMolFile());
		}
		MolAtom attachPoint = mol.getAtom(atomId - 1);
		if (attachPoint.getAtno() == 1) {
			throw new Exception("Hydrogen can't be used as attaching point.");
		}
		int explictHCount = 0;
		for (int i = 0; i < attachPoint.getBondCount(); i++) {
			MolBond bond = attachPoint.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(attachPoint);
			if (neighborAtom.getAtno() == 1) {
				explictHCount += 1;
			}
		}

		if ((attachPoint.getImplicitHcount() == 0 && explictHCount == 0) || attachPoint.getBondCount() == 0) {
			throw new Exception("No free valence to add new group.");
		}
		MolAtom removedHatom = null;
		for (int i = 0; i < attachPoint.getBondCount(); i++) {
			MolBond bond = attachPoint.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(attachPoint);
			if (neighborAtom.getAtno() == 1) {
				mol.removeBond(bond);
				mol.removeAtom(neighborAtom);
				removedHatom = neighborAtom;
				break;
			}
		}

		OEUnaryAtomPred pred = null;
		Vector<Integer> idList = new Vector<Integer>();
		for (MolAtom atom : mol.getAtomArray()) {
			if (atom.getAtno() > 1) {
				idList.add(mol.indexOf(atom));
			}
		}
		MolAtom attachPoint2 = null;
		for (MolAtom atom : frag.getAtomArray()) {
			if (frag.indexOf(atom) == 0) {
				attachPoint2 = atom;
			}
			mol.insertAtom(mol.getAtomCount(), atom);

		}

		if (attachPoint2 == null) {
			throw new Exception("No Attach point defined in the group.");
		}

		if (removedHatom != null) {
			CTransform3D transform3D = new CTransform3D();
			transform3D.setTranslation(removedHatom.getX(), removedHatom.getY(), removedHatom.getZ());
			frag.transform(transform3D);
		}

		for (MolBond bond : frag.getBondArray()) {
			mol.insertBond(mol.getBondCount(), bond);
		}

		MolBond bond = new MolBond(attachPoint, attachPoint2);

		mol.insertBond(mol.getBondCount(), bond);

		OEGraphMol oemol = OEChemFunc.getInstance().convertChemAxonMol(mol);
		oechem.OEAddExplicitHydrogens(oemol);
		Vector<Integer> ids = new Vector<>();
		for (OEAtomBaseIter iter = oemol.GetAtoms(); iter.hasNext(); ) {
			OEAtomBase atom = iter.next();
			if (idList.contains(atom.GetIdx())) {
				Integer idx = atom.GetIdx();
				ids.add(idx);
//				if (pred == null) {
//					pred = new OEHasAtomIdx(idx);
//				} else {
//					pred = new OEOrAtom(pred, new OEHasAtomIdx(idx));
//				}
			}
		}
		oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol);
		Protonator.getInstance().protonate(oemol);
		OEGraphMol minimizedMol = OEChemFunc.getInstance().minimizeOEMol(oemol, ids);

		return OEChemFunc.getInstance().convertOEChemMol(minimizedMol).toFormat("mol");
	}

	public static String deleteAtomFromMolecule(String molString, int atomId) throws Exception {
		Molecule mol = MolImporter.importMol(molString);
		MolAtom atomToBeDeleted = mol.getAtom(atomId - 1);
		for (int i = 0; i < atomToBeDeleted.getBondCount(); i++) {
			MolBond bond = atomToBeDeleted.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(atomToBeDeleted);
			mol.removeBond(bond);
			if (neighborAtom.getAtno() == 1) {
				mol.removeAtom(neighborAtom);
			}
		}
		mol.removeAtom(atomToBeDeleted);

		return OEChemFunc.getInstance().fixHydrogen(mol).toFormat("mol");
	}

	public static String linkAtoms(String molString, int atomId1, int atomId2, int length) throws Exception {

		Molecule mol = MolImporter.importMol(molString);
		MolAtom headAtom = mol.getAtom(atomId1 - 1);
		MolAtom tailAtom = mol.getAtom(atomId2 - 1);
		if (headAtom == null || tailAtom == null) {
			throw new Exception("No Atom found.");
		}


		MolAtom attachPoint1 = mol.getAtom(atomId1 - 1);
		if (attachPoint1.getAtno() == 1) {
			throw new Exception("Hydrogen can't be used as attaching point.");
		}

		int explictHCount = 0;
		for (int i = 0; i < attachPoint1.getBondCount(); i++) {
			MolBond bond = attachPoint1.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(attachPoint1);
			if (neighborAtom == tailAtom && length == 0) {
				throw new Exception("Two atoms are already linked.");
			}
			if (neighborAtom.getAtno() == 1) {
				explictHCount += 1;
			}
		}

		if ((attachPoint1.getImplicitHcount() == 0 && explictHCount == 0) || attachPoint1.getBondCount() == 0) {
			throw new Exception("No free valence to add new group.");
		}

		MolAtom attachPoint2 = mol.getAtom(atomId2 - 1);
		if (attachPoint2.getAtno() == 1) {
			throw new Exception("Hydrogen can't be used as attaching point.");
		}

		explictHCount = 0;
		for (int i = 0; i < attachPoint2.getBondCount(); i++) {
			MolBond bond = attachPoint2.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(attachPoint2);
			if (neighborAtom.getAtno() == 1) {
				explictHCount += 1;
			}
		}

		if ((attachPoint2.getImplicitHcount() == 0 && explictHCount == 0) || attachPoint2.getBondCount() == 0) {
			throw new Exception("No free valence to add new group.");
		}


		MolAtom removedHatom1 = null;
		for (int i = 0; i < attachPoint1.getBondCount(); i++) {
			MolBond bond = attachPoint1.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(attachPoint1);
			if (neighborAtom.getAtno() == 1) {
				mol.removeBond(bond);
				mol.removeAtom(neighborAtom);
				removedHatom1 = neighborAtom;
				break;
			}
		}

		MolAtom removedHatom2 = null;
		for (int i = 0; i < attachPoint2.getBondCount(); i++) {
			MolBond bond = attachPoint2.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(attachPoint2);
			if (neighborAtom.getAtno() == 1) {
				mol.removeBond(bond);
				mol.removeAtom(neighborAtom);
				removedHatom2 = neighborAtom;
				break;
			}
		}

		if (length == -1) {
			length = getPreferredLength(attachPoint1, attachPoint2);
		}

		if (length > 0) {
			OEUnaryAtomPred pred = null;
			Vector<Integer> idList = new Vector<Integer>();
			for (MolAtom atom : mol.getAtomArray()) {
				if (atom.getAtno() > 1) {
					idList.add(mol.indexOf(atom));
				}
			}
			MolAtom attachPoint3 = null;
			MolAtom attachPoint4 = null;
			StringBuilder sb = new StringBuilder();
			for(int i=0;i<length;i++){
				sb.append("C");
			}
			String fragSmi = sb.toString();

			//todo: write a function to get 3d conformation of the linker.
			Molecule fragMol = MolImporter.importMol(fragSmi);
			fragMol.clean(3, null);
			//Molecule fragMol = ModelingTasks.getCorina3DStructure(fragSmi).get(0);
			for (MolAtom atom : fragMol.getAtomArray()) {
				if (fragMol.indexOf(atom) == 0) {
					attachPoint3 = atom;
				}
				if (fragMol.indexOf(atom) == length - 1) {
					attachPoint4 = atom;
				}
				mol.insertAtom(mol.getAtomCount(), atom);
			}

			if (attachPoint3 == null || attachPoint4 == null) {
				throw new Exception("No Attach point defined in the group.");
			}

			if (removedHatom1 != null && removedHatom2 != null) {
				CTransform3D transform3D = new CTransform3D();
				transform3D.setTranslation((removedHatom1.getX() + removedHatom2.getX()) / 2.0, (removedHatom1.getY() + removedHatom2.getY()) / 2.0, (removedHatom1.getZ() + removedHatom2.getZ()) / 2.0);
				fragMol.transform(transform3D);
			}


			if (removedHatom1 != null) {
				attachPoint3.setXYZ(removedHatom1.getX(), removedHatom1.getY(), removedHatom1.getZ());
			}

			if (removedHatom2 != null) {
				attachPoint4.setXYZ(removedHatom2.getX(), removedHatom2.getY(), removedHatom2.getZ());
			}


			for (MolBond bond : fragMol.getBondArray()) {
				mol.insertBond(mol.getBondCount(), bond);
			}

			MolBond bond1 = new MolBond(attachPoint1, attachPoint3);
			mol.insertBond(mol.getBondCount(), bond1);

			MolBond bond2 = new MolBond(attachPoint2, attachPoint4);
			mol.insertBond(mol.getBondCount(), bond2);

			OEGraphMol oemol = OEChemFunc.getInstance().convertChemAxonMol(mol);
			Vector<Integer> ids = new Vector<>();
			for (OEAtomBaseIter iter = oemol.GetAtoms(); iter.hasNext(); ) {
				OEAtomBase atom = iter.next();
				if (idList.contains(atom.GetIdx())) {
					Integer idx = atom.GetIdx();
					ids.add(idx);
//					if (pred == null) {
//						pred = new OEHasAtomIdx(idx);
//					} else {
//						pred = new OEOrAtom(pred, new OEHasAtomIdx(idx));
//					}
				}
			}
			oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol);
			oechem.OESuppressHydrogens(oemol);
			oechem.OEAddExplicitHydrogens(oemol);
			Protonator.getInstance().protonate(oemol);
			OEGraphMol minimizedMol = OEChemFunc.getInstance().minimizeOEMol(oemol, ids);
			Molecule molecule = OEChemFunc.getInstance().convertOEChemMol(minimizedMol);
			return molecule.toFormat("mol");
			//return mol.toFormat("mol");
		} else {
			MolBond bond = new MolBond(attachPoint1, attachPoint2);
			mol.insertBond(mol.getBondCount(), bond);
			return mol.toFormat("mol");
		}


	}

	public static String modifyAtom(String molString, int atomId, MorphingElement element) throws Exception {
		Molecule mol = MolImporter.importMol(molString);
		MolAtom atomToBeChanged = mol.getAtom(atomId - 1);
		atomToBeChanged.setAtno(element.getElementNo());
		for (int i = 0; i < atomToBeChanged.getBondCount(); i++) {
			MolBond bond = atomToBeChanged.getBond(i);
			MolAtom neighborAtom = bond.getOtherAtom(atomToBeChanged);
			if (neighborAtom.getAtno() == 1) {
				mol.removeBond(bond);
				mol.removeAtom(neighborAtom);
			}
		}
		return OEChemFunc.getInstance().fixHydrogen(mol).toFormat("mol");
	}


	public static String modifyBond(String molString, int idx1, int idx2, String morphing_bond_type) throws Exception {
		Molecule mol = MolImporter.importMol(molString);
		MolAtom headAtom = mol.getAtom(idx1 - 1);
		MolAtom tailAtom = mol.getAtom(idx2 - 1);
		if (headAtom == null || tailAtom == null) {
			throw new Exception("No Atom found.");
		}
		MolBond bondToModify = null;
		for (MolBond bond : mol.getBondArray()) {
			if ((bond.getAtom1() == headAtom && bond.getAtom2() == tailAtom) || (bond.getAtom1() == tailAtom && bond.getAtom2() == headAtom)) {
				bondToModify = bond;
				break;
			}
		}
		if (bondToModify == null) {
			throw new Exception("No Bond found.");
		}

		int type = 1;
		if (morphing_bond_type.equals("single")) {
			type = 1;
			bondToModify.setType(type);
		} else if (morphing_bond_type.equals("double")) {
			type = 2;
			bondToModify.setType(type);
		} else if (morphing_bond_type.equals("triple")) {
			type = 3;
			bondToModify.setType(type);
		} else if (morphing_bond_type.equals("none")) {
			mol.removeBond(bondToModify);
		} else {
			throw new Exception("Unknown bond type.");
		}

		return OEChemFunc.getInstance().fixHydrogen(mol).toFormat("mol");

	}
}
