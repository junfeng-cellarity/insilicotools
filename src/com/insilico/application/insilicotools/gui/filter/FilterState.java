package com.insilico.application.insilicotools.gui.filter;

public abstract class FilterState {
    private boolean isNewFilter = true;

    public boolean isNewFilter() {
        return isNewFilter;
    }

    public void setNewFilter(boolean newFilter) {
        isNewFilter = newFilter;
    }

    public abstract String isValidState();
}