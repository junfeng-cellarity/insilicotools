package com.insilico.application.insilicotools.gui.filter.fragment;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

/**
 * Created by IntelliJ IDEA.
 * User: liuh
 * Date: 6/26/12
 * Time: 4:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class FragmentFilterTransport extends FilterTransport {

    public FragmentFilterTransport(final FilterState state) {
        super(state);
    }

    public FilterController getFilterController(final String name, final TreeFilter treeFilter) {
        return new FragmentFilterController(name, state, treeFilter);
    }
}
