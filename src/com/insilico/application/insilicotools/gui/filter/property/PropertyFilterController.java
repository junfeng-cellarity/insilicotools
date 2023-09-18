package com.insilico.application.insilicotools.gui.filter.property;

import com.insilico.application.insilicotools.gui.filter.*;

public class PropertyFilterController extends GuiFilterController {
    public PropertyFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new PropertyFilterState(), treeFilter);
    }

    public PropertyFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    public FilterTransport getTransport() {
        return new PropertyFilterTransport(state);
    }

    protected Filter getFilter() {
        return new PropertyFilter();
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new PropertyFilterGui((PropertyFilterState) state);
    }
}