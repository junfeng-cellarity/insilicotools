package com.insilico.application.insilicotools.gui.filter.fragment;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class FragmentFilterFactory extends FilterFactory {
    public FragmentFilterFactory() {
        super("Fragment");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new FragmentFilterController(getName(), treeFilter);
    }
}