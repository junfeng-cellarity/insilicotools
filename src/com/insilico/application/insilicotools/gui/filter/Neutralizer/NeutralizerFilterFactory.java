package com.insilico.application.insilicotools.gui.filter.Neutralizer;

import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class NeutralizerFilterFactory extends FilterFactory {
    public NeutralizerFilterFactory() {
        super("Neutralizer");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new NeutralizerFilterController(getName(),treeFilter);
    }
}
