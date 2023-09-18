package com.insilico.application.insilicotools.gui.filter.diverse;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class DiverseFilterFactory extends FilterFactory {
    public DiverseFilterFactory() {
        super("Diverse");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new DiverseFilterController(getName(), treeFilter);
    }
}
