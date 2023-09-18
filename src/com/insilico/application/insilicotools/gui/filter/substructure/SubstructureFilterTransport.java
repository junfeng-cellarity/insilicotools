package com.insilico.application.insilicotools.gui.filter.substructure;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class SubstructureFilterTransport extends FilterTransport {
    public SubstructureFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new SubstructureFilterController(name, state, treeFilter);
    }
}
