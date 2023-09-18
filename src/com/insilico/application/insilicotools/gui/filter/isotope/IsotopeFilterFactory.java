package com.insilico.application.insilicotools.gui.filter.isotope;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class IsotopeFilterFactory extends FilterFactory {
    public IsotopeFilterFactory() {
        super("Isotope");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new IsotopeFilterController(getName(), treeFilter);
    }
}