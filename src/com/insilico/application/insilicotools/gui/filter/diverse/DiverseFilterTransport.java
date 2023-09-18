package com.insilico.application.insilicotools.gui.filter.diverse;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class DiverseFilterTransport extends FilterTransport {
    public DiverseFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new DiverseFilterController(name, state, treeFilter);
    }
}
