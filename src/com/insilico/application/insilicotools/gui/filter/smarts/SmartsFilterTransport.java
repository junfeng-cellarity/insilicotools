package com.insilico.application.insilicotools.gui.filter.smarts;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class SmartsFilterTransport extends FilterTransport {
    public SmartsFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new SmartsFilterController(name, state, treeFilter);
    }
}
