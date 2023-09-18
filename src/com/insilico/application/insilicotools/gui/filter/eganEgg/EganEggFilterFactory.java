package com.insilico.application.insilicotools.gui.filter.eganEgg;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class EganEggFilterFactory extends FilterFactory {
    public EganEggFilterFactory() {
        super("Egan Egg");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new EganEggFilterController(getName(), treeFilter);
    }
}