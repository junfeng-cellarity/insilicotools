package com.insilico.application.insilicotools.gui.filter.substructure;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class SubstructureFilterFactory extends FilterFactory {
    public SubstructureFilterFactory() {
        super("Substructure");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new SubstructureFilterController(getName(), treeFilter);
    }
}