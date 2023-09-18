package com.insilico.application.insilicotools.gui.filter;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;

public abstract class Filter {

    public abstract FilterResult filter(final ProgressReporter reporter, final PropertyMolecule[] mols,
                                        final FilterState state) throws Exception;

}