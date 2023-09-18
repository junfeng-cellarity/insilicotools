package com.insilico.application.insilicotools.gui.filter.Neutralizer;

import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;
import com.insilico.application.insilicotools.gui.filter.FilterState;

public class NeutralizerFilterTransport extends FilterTransport {
    public NeutralizerFilterTransport(FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new NeutralizerFilterController(name,state,treeFilter);
    }
}
