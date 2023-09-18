package com.insilico.application.insilicotools.gui.filter.eganEgg;

import com.insilico.application.insilicotools.gui.filter.*;

public class EganEggFilterController extends GuiFilterController {
    public EganEggFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new EganEggFilterState(), treeFilter);
    }

    public EganEggFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected Filter getFilter() {
        return new EganEggFilter();
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new EganEggFilterGui((EganEggFilterState) state);
    }

    public FilterTransport getTransport() {
        return new EganEggFilterTransport(state);
    }
}