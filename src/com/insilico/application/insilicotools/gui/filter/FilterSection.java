package com.insilico.application.insilicotools.gui.filter;


import com.insilico.application.insilicotools.data.PropertyMolecule;

public class FilterSection {
    private final PropertyMolecule[] molecules;
    private final String text;
    private final String tooltip;
    private final int state;

    public final static int PASS = 0;
    public final static int FAIL = 1;
    public final static int NEITHER = 2;

    public FilterSection(final PropertyMolecule[] molecules, final String text) {
        this(molecules, text, null, NEITHER);
    }

    public FilterSection(final PropertyMolecule[] molecules, final String text, final int state) {
        this(molecules, text, null, state);
    }

    public FilterSection(final PropertyMolecule[] molecules, final String text, final String tooltip, final int state) {
        this.molecules = molecules;
        this.text = text;
        this.tooltip = tooltip;
        this.state = state;
    }

    public PropertyMolecule[] getMolecules() {
        return molecules;
    }

    public String getText() {
        return text;
    }

    public String getTooltip() {
        return tooltip;
    }

    public int getState() {
        return state;
    }

    public String toString() {
        return text;
    }
}