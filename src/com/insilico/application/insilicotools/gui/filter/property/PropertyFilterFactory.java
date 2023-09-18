package com.insilico.application.insilicotools.gui.filter.property;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class PropertyFilterFactory extends FilterFactory {
    public PropertyFilterFactory() {
        super("Properties");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new PropertyFilterController(getName(), treeFilter);
    }
}