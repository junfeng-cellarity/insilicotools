package com.insilico.application.insilicotools.gui.filter.uniq;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class UniqFilterFactory extends FilterFactory{
    public UniqFilterFactory() {
        super("Unique");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
     return new UniqFilterController(getName(),treeFilter);
    }
}
