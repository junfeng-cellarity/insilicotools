package com.insilico.application.insilicotools.util;

import chemaxon.struc.Molecule;
import chemaxon.struc.DPoint3;
import chemaxon.struc.CTransform3D;
import chemaxon.struc.MolAtom;
import chemaxon.util.MolHandler;

import java.util.Vector;

public class RigidBodySuperposition {
    private chemaxon.struc.Molecule molTemplate;
    private chemaxon.struc.Molecule molTarget;
    private Vector idTemplate;
    private Vector idTarget;
    private int [] alignArray = null;
    double totalRMSD = -1;

    public RigidBodySuperposition(chemaxon.struc.Molecule molTemplate, chemaxon.struc.Molecule molTarget, Vector idTemplate, Vector idTarget) {
        this.molTemplate = molTemplate;
        this.molTarget = molTarget;
        this.idTemplate = idTemplate;
        this.idTarget = idTarget;
    }

    public RigidBodySuperposition(Molecule molTemplate, Molecule molTarget, int [] match){
        this.molTemplate = molTemplate;
        this.molTarget = molTarget.cloneMolecule();
        //if(match.length!=molTemplate.getAtomCount()) return;
        idTemplate = new Vector();
        idTarget = new Vector();
        for (int i = 0; i < match.length; i++) {
            int id = match[i];
            if(id==-1) continue;
            idTemplate.add(new Integer(i+1));
            idTarget.add(new Integer(match[i]+1));
        }
        alignArray = match;

    }

    public Molecule superimpose() {
        if(idTemplate.size()!=idTarget.size()){
            System.out.println("Superposition requires target and template are the same size.");
            return null;
        }

        if(idTemplate.size()<3){
            System.out.println("Superposition requires at least 3 pairs of atoms specified.");
            return null;
        }

        int n = molTemplate.getAtomCount();

        if(alignArray==null){
            alignArray = new int [n];
            for (int i = 0; i < alignArray.length; i++) {
                if(idTemplate.contains(new Integer(i+1))){
                    alignArray[i] = ((Integer)idTarget.get(idTemplate.indexOf(new Integer(i+1)))).intValue()-1;
                }else{
                    alignArray[i] = -1;
                }
            }
        }

        Molecule newMol = SuperImpose(alignArray);

/*        if(newMol.getAtomCount()==molTemplate.getAtomCount()){
            totalRMSD = 0.0;
            for(int i=0;i<newMol.getAtomCount();i++){
                MolAtom at1 = newMol.getAtom(i);
                MolAtom at2 = molTemplate.getAtom(i);
                totalRMSD += (at1.getX()-at2.getX())*(at1.getX()-at2.getX())+(at1.getY()-at2.getY())*(at1.getY()-at2.getY())+(at1.getZ()-at2.getZ())*(at1.getZ()-at2.getZ());
            }
            totalRMSD = Math.sqrt(totalRMSD/newMol.getAtomCount());
        }else{ */
        totalRMSD = 0.0;
        int count = 0;
        for(int i=0;i<alignArray.length;i++){
            if(alignArray[i]==-1){
                continue;
            }else{
                MolAtom at1 = molTemplate.getAtom(i);
                MolAtom at2 = newMol.getAtom(alignArray[i]);
                totalRMSD += (at1.getX()-at2.getX())*(at1.getX()-at2.getX())+(at1.getY()-at2.getY())*(at1.getY()-at2.getY())+(at1.getZ()-at2.getZ())*(at1.getZ()-at2.getZ());
                count += 1;
            }
        }
        totalRMSD = Math.sqrt(totalRMSD/count);
//        }
        return newMol;

    }

    public double getTotalRMSD() {
        return totalRMSD;
    }

    private Molecule SuperImpose(int[] alignArray) {
        MolHandler mh = new MolHandler(molTemplate);
        mh.align(molTarget,alignArray);

        DPoint3 centerTarget = calc_center(molTarget,idTarget);
        DPoint3 centerTemplate = calc_center(molTemplate,idTemplate);

        double deltax = centerTemplate.x - centerTarget.x;
        double deltay = centerTemplate.y - centerTarget.y;
        double deltaz = centerTemplate.z - centerTarget.z;

        for(int i=0;i<molTarget.getAtomArray().length;i++){
            double x = molTarget.getAtomArray()[i].getX();
            double y = molTarget.getAtomArray()[i].getY();
            double z = molTarget.getAtomArray()[i].getZ();
            molTarget.getAtomArray()[i].setX(x+deltax);
            molTarget.getAtomArray()[i].setY(y+deltay);
            molTarget.getAtomArray()[i].setZ(z+deltaz);
        }

//        CTransform3D translation = new CTransform3D();
//        translation.setTranslation(centerTemplate.x-centerTarget.x,centerTemplate.y-centerTarget.y,centerTemplate.z-centerTarget.z);
//
//        molTarget.transform(translation);

        return molTarget;
    }

    private DPoint3 calc_center(Molecule mol, Vector idList){
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        for(int i =0;i<mol.getAtomCount();i++){
            Integer id = new Integer(i+1);
            if(idList.contains(id)){
                x += mol.getAtom(i).getX();
                y += mol.getAtom(i).getY();
                z += mol.getAtom(i).getZ();
            }
        }
        x/=(double)idList.size();
        y/=(double)idList.size();
        z/=(double)idList.size();

        return new DPoint3(x,y,z);

    }

}
