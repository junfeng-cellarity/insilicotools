package com.insilico.application.insilicotools.gui.filter.uniq;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class UniqFilterTransport extends FilterTransport{
    protected UniqFilterTransport(FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new UniqFilterController(name,state,treeFilter);  //To change body of implemented methods use File | Settings | File Templates.
    }
}
