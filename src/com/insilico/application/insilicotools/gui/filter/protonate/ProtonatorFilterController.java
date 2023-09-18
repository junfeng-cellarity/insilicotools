package com.insilico.application.insilicotools.gui.filter.protonate;

import com.insilico.application.insilicotools.gui.filter.*;

public class ProtonatorFilterController extends GuiFilterController {
    public ProtonatorFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new ProtonatorFilterState(), treeFilter);
    }

    protected ProtonatorFilterController(String name, FilterState state, TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new ProtonatorFilterGui((ProtonatorFilterState)state);
    }

    protected Filter getFilter() {
        return new ProtonatorFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    public FilterTransport getTransport(){
        return new ProtonatorFilterTransport(state);
    }

}
