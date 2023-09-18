package com.insilico.application.insilicotools.gui.filter.names;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class NameFilterTransport extends FilterTransport {
    public NameFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new NameFilterController(name, state, treeFilter);
    }
}
