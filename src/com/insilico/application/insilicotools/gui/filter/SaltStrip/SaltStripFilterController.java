package com.insilico.application.insilicotools.gui.filter.SaltStrip;

import com.insilico.application.insilicotools.gui.filter.*;

public class SaltStripFilterController extends GuiFilterController {
    public SaltStripFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new SaltStripFilterState(), treeFilter);
    }

    protected SaltStripFilterController(String name, FilterState state, TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new SaltStripFilterGui((SaltStripFilterState)state);  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected Filter getFilter() {
        return new SaltStripFilter();  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    public FilterTransport getTransport(){
          return new SaltStripFilterTransport(state);
    }
}
