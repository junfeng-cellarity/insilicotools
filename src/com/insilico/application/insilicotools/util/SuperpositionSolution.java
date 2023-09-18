package com.insilico.application.insilicotools.util;

import openeye.oechem.OEGraphMol;

public class SuperpositionSolution implements Comparable {
    double strainEnergy;
    double shapeScore;

    OEGraphMol template;
    OEGraphMol target;

    public SuperpositionSolution(OEGraphMol template, OEGraphMol target, double strainEnergy, double shapeScore) {
        this.template = template;
        this.target = target;
        this.strainEnergy = strainEnergy;
        this.shapeScore = shapeScore;
    }

    public OEGraphMol getTemplate() {
        return template;
    }

    public OEGraphMol getTarget() {
        return target;
    }

    public int compareTo(Object o) {
        Double score = this.getScore();
        Double score2 = ((SuperpositionSolution) o).getScore();
        return score.compareTo(score2);
    }


    public double getStrainEnergy() {
        return strainEnergy;
    }

    public void setStrainEnergy(double strainEnergy) {
        this.strainEnergy = strainEnergy;
    }

    public double getShapeScore() {
        return shapeScore;
    }

    public void setShapeScore(double shapeScore) {
        this.shapeScore = shapeScore;
    }

    public double getScore() {
        return shapeScore;
    }

    public double getEnergy() {
        return strainEnergy;
    }

    public String getXAxisName() {
        return "Strain Energy";
    }

    public String getYAxisName() {
        return "Shape/Feature Score";
    }

}
