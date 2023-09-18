package com.insilico.application.insilicotools.gui.filter.Neutralizer;

import com.insilico.application.insilicotools.gui.filter.*;

public class NeutralizerFilterController extends GuiFilterController {
    public NeutralizerFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new NeutralizerFilterState(), treeFilter);
    }

    protected NeutralizerFilterController(String name, FilterState state, TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected FilterGui getGuiComponent(FilterState state) {
        return new NeutralizerFilterGui((NeutralizerFilterState)state);  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected Filter getFilter() {
        return new NeutralizerFilter();  //To change body of implemented methods use File | Settings | File Templates.
    }

    protected boolean canUseDeprotectedSmiles() {
        return false;
    }

    public FilterTransport getTransport(){
          return new NeutralizerFilterTransport(state);
    }
}
