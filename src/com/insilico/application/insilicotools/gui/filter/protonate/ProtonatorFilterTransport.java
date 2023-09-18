package com.insilico.application.insilicotools.gui.filter.protonate;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.Neutralizer.NeutralizerFilterController;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class ProtonatorFilterTransport extends FilterTransport {
    public ProtonatorFilterTransport(FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new ProtonatorFilterController(name,state,treeFilter);
    }
}
