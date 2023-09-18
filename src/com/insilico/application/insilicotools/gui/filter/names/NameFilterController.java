package com.insilico.application.insilicotools.gui.filter.names;

import com.insilico.application.insilicotools.gui.filter.*;

public class NameFilterController extends GuiFilterController {
    public NameFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new NameFilterState(), treeFilter);
    }

    public NameFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected Filter getFilter() {
        return new NameFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new NameFilterGui((NameFilterState) state);
    }

    public FilterTransport getTransport() {
        return new NameFilterTransport(state);
    }
}