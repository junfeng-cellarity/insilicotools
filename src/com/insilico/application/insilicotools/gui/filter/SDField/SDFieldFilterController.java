package com.insilico.application.insilicotools.gui.filter.SDField;

import com.insilico.application.insilicotools.gui.filter.*;

import java.util.Vector;

public class SDFieldFilterController extends GuiFilterController {
    Vector<String> properties;

    public SDFieldFilterController(final String name, final TreeFilter treeFilter, Vector<String> properties) {
        this(name, new SDFieldFilterState(), treeFilter);
        this.properties = properties;
    }

    public SDFieldFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
        this.properties = ((SDFieldFilterState)state).getPropertyKeys();
    }

    public FilterTransport getTransport() {
        return new SDFieldFilterTransport(state);
    }

    protected Filter getFilter() {
        return new SDFieldFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    protected FilterGui getGuiComponent(FilterState state) {
        SDFieldFilterState filterState = (SDFieldFilterState) state;
        filterState.setPropertyKeys(properties);
        return new SDFieldFilterGui(filterState);
    }
}
