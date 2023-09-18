package com.insilico.application.insilicotools.gui.filter;

import java.util.EventListener;

public interface MoleculeSelectionListener extends EventListener {
    public void moleculesSelected(MoleculeSelectionEvent evt);
}