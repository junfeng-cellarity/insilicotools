package com.insilico.application.insilicotools.gui.filter;

import com.insilico.application.insilicotools.data.PropertyMolecule;

import java.util.EventObject;

public class MoleculeSelectionEvent extends EventObject {
    private final PropertyMolecule[] mols;
    private final boolean isValueAdjusting;

    public MoleculeSelectionEvent(Object source, PropertyMolecule[] mols, boolean isValueAdjusting) {
        super(source);
        this.mols = mols;
        this.isValueAdjusting = isValueAdjusting;
    }

    public PropertyMolecule[] getMols() {
        return mols;
    }

    public boolean isValueAdjusting() {
        return isValueAdjusting;
    }
}