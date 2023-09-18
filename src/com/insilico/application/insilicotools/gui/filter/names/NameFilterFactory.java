package com.insilico.application.insilicotools.gui.filter.names;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class NameFilterFactory extends FilterFactory {
    public NameFilterFactory() {
        super("Names");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new NameFilterController(getName(), treeFilter);
    }
}