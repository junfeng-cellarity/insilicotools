package com.insilico.application.insilicotools.gui.filter.smarts;

import com.insilico.application.insilicotools.gui.filter.*;

public class SmartsFilterController extends GuiFilterController {
    public SmartsFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new SmartsFilterState(), treeFilter);
    }

    public SmartsFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected Filter getFilter() {
        return new SmartsFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new SmartsFilterGui((SmartsFilterState) state);
    }

    public FilterTransport getTransport() {
        return new SmartsFilterTransport(state);
    }
}