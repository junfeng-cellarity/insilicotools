package com.insilico.application.insilicotools.gui.filter.eganEgg;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class EganEggFilterTransport extends FilterTransport {
    public EganEggFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new EganEggFilterController(name, state, treeFilter);
    }
}
