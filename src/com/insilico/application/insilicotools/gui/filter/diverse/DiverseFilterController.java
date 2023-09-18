package com.insilico.application.insilicotools.gui.filter.diverse;

import com.insilico.application.insilicotools.gui.filter.*;

public class DiverseFilterController extends GuiFilterController {
    public DiverseFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new DiverseFilterState(), treeFilter);
    }

    public DiverseFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected Filter getFilter() {
        return new DiverseFilter();
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new DiverseFilterGui((DiverseFilterState) state);
    }

    public FilterTransport getTransport() {
        return new DiverseFilterTransport(state);
    }
}