package com.insilico.application.insilicotools.gui.filter.smarts;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class SmartsFilterFactory extends FilterFactory {
    public SmartsFilterFactory() {
        super("Smarts");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new SmartsFilterController(getName(), treeFilter);
    }
}
