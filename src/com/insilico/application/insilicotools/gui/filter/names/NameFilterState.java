package com.insilico.application.insilicotools.gui.filter.names;

import com.insilico.application.insilicotools.gui.filter.FilterState;

public class NameFilterState extends FilterState {
    private String names;

    public String isValidState() {
        return names == null || names.trim().length() == 0 ? "You must enter a list of names" : null;
    }

    public String getNames() {
        return names;
    }

    public void setNames(final String names) {
        this.names = names;
    }
}