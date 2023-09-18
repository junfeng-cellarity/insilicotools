package com.insilico.application.insilicotools.gui.filter.isotope;

import com.insilico.application.insilicotools.gui.filter.*;

public class IsotopeFilterController extends FilterController {
    public IsotopeFilterController(final String name, final TreeFilter treeFilter) {
        super(name, null, treeFilter);
    }
    public IsotopeFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
            super(name, state, treeFilter);
        }

    protected Filter getFilter() {
        return new IsotopeFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    public FilterTransport getTransport() {
        return new IsotopeFilterTransport(state);
    }
}