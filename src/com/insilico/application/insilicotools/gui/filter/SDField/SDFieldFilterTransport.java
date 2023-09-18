package com.insilico.application.insilicotools.gui.filter.SDField;

import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;
import com.insilico.application.insilicotools.gui.filter.FilterState;

public class SDFieldFilterTransport extends FilterTransport {
   public SDFieldFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new SDFieldFilterController(name, state, treeFilter);
    }
}
