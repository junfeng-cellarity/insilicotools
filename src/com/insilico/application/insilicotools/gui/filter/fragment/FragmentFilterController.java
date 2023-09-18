package com.insilico.application.insilicotools.gui.filter.fragment;

import com.insilico.application.insilicotools.gui.filter.*;


public class FragmentFilterController extends FilterController {
    public FragmentFilterController(final String name, final TreeFilter treeFilter) {
        super(name, null, treeFilter);
    }
    public FragmentFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected Filter getFilter() {
        return new FragmentFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }
     public FilterTransport getTransport() {
        return new FragmentFilterTransport(state);
    }
}