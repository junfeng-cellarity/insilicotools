package com.insilico.application.insilicotools.gui.filter.uniq;

import com.insilico.application.insilicotools.gui.filter.*;

public class UniqFilterController extends GuiFilterController {
    public UniqFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new UniqFilterState(), treeFilter);
    }

    protected UniqFilterController(String name, FilterState state, TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new UniqFilterGui((UniqFilterState)state);  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected Filter getFilter() {
        return new UniqFilter();  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected boolean canUseDeprotectedSmiles() {
        return true;
    }

    public FilterTransport getTransport(){
          return new UniqFilterTransport(state);
    }
}
