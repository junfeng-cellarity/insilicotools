package com.insilico.application.insilicotools.gui.filter.isotope;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

/**
 * Created by IntelliJ IDEA.
 * User: liuh
 * Date: 6/28/12
 * Time: 12:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class IsotopeFilterTransport extends FilterTransport {
    protected IsotopeFilterTransport(FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new IsotopeFilterController(name,state,treeFilter);
    }
}
