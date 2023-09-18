package com.insilico.application.insilicotools.gui.filter.alert;

import com.insilico.application.insilicotools.gui.filter.*;

public class AlertFilterController extends GuiFilterController {
    public AlertFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new AlertFilterState(), treeFilter);
    }

    public AlertFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    public FilterTransport getTransport() {
        return new AlertFilterTransport(state);
    }

    protected Filter getFilter() {
        return new AlertFilter();
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new AlertFilterGui((AlertFilterState) state);
    }
}