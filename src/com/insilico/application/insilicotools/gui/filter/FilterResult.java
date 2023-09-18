package com.insilico.application.insilicotools.gui.filter;

import com.insilico.application.insilicotools.data.PropertyMolecule;

import java.util.ArrayList;
import java.util.Arrays;

public class FilterResult {
    private final ArrayList<FilterSection> sections = new ArrayList<FilterSection>();

    public FilterResult(final PropertyMolecule[] pass, final PropertyMolecule[] fail) {
        if(pass != null) {
            sections.add(new FilterSection(pass, "Pass"));
        }
        else {
            throw new NullPointerException("Passed molecule array cannot be null");
        }
        if(fail != null) {
            sections.add(new FilterSection(fail, "Fail"));
        }
        else {
            throw new NullPointerException("Fail molecule array cannot be null");
        }
    }

    public FilterResult(FilterSection[] sections) {
        this.sections.addAll(Arrays.asList(sections));
    }

    public FilterSection[] getSections() {
        return sections.toArray(new FilterSection[sections.size()]);
    }
}