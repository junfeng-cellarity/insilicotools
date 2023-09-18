package com.insilico.application.insilicotools.gui.filter.SaltStrip;

import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class SaltStripFilterFactory extends FilterFactory {
    public SaltStripFilterFactory() {
        super("Salt Strip");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
     return new SaltStripFilterController(getName(),treeFilter);
    }
}
