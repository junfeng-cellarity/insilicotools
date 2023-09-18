package com.insilico.application.insilicotools.gui.filter.SaltStrip;

import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class SaltStripFilterTransport extends FilterTransport {
    protected SaltStripFilterTransport(FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new SaltStripFilterController(name,state,treeFilter);
    }
}
