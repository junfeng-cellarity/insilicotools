package com.insilico.application.insilicotools.gui.filter.protonate;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.Neutralizer.NeutralizerFilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class ProtonatorFilterFactory extends FilterFactory {
    public ProtonatorFilterFactory() {
        super("Protonator");
    }

    public FilterController getInstance(final TreeFilter treeFilter) {
        return new ProtonatorFilterController(getName(),treeFilter);
    }

}
