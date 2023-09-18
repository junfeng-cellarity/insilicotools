package com.insilico.application.insilicotools.gui.filter;

public abstract class FilterFactory {
    private String name;

    public FilterFactory() {
    }

    protected FilterFactory(final String name) {
        this.name = name;
    }

    public abstract FilterController getInstance(final TreeFilter treeFilter);

    public String getName() {
        return name;
    }

    public void setName(String name){
        this.name = name;
    }

    public String toString() {
        return getName();
    }
}