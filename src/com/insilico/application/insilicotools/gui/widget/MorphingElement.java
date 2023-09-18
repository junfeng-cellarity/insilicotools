package com.insilico.application.insilicotools.gui.widget;

import java.util.Vector;

public class MorphingElement {
    int elementNo;
    String elementSymbol;

    public MorphingElement(int elementNo, String elementSymbol) {
        this.elementNo = elementNo;
        this.elementSymbol = elementSymbol;

    }

    public int getElementNo() {
        return elementNo;
    }

    public void setElementNo(int elementNo) {
        this.elementNo = elementNo;
    }

    public String getElementSymbol() {
        return elementSymbol;
    }

    public void setElementSymbol(String elementSymbol) {
        this.elementSymbol = elementSymbol;
    }

    @Override
    public String toString() {
        return elementSymbol;
    }

    public static Vector<MorphingElement> getCommonElements() {
        Vector<MorphingElement> elements = new Vector<MorphingElement>();
        elements.add(new MorphingElement(1, "H"));
        elements.add(new MorphingElement(6, "C"));
        elements.add(new MorphingElement(7, "N"));
        elements.add(new MorphingElement(8, "O"));
        elements.add(new MorphingElement(9, "F"));
        elements.add(new MorphingElement(15, "P"));
        elements.add(new MorphingElement(16, "S"));
        elements.add(new MorphingElement(17, "Cl"));
        elements.add(new MorphingElement(35, "Br"));
        elements.add(new MorphingElement(53, "I"));

        return elements;
    }
}
