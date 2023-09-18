package com.insilico.application.insilicotools.gui.filter;

public abstract class FilterTransport {
    protected final FilterState state;

    protected FilterTransport(final FilterState state) {
        this.state = state;
    }

    protected FilterState getState() {
        return state;
    }

    public abstract FilterController getFilterController(final String name, final TreeFilter treeFilter);

}