package com.insilico.application.insilicotools.gui.filter.substructure;

import com.insilico.application.insilicotools.gui.filter.FilterState;

public class SubstructureFilterState extends FilterState {
    private String query;

    public String getQuery() {
        return query;
    }

    public void setQuery(final String query) {
        this.query = query;
    }

    public String isValidState() {
        return query == null ? "A substructure must be defined" : null;
    }
}