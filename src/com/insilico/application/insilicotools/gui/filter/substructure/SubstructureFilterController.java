package com.insilico.application.insilicotools.gui.filter.substructure;

import com.insilico.application.insilicotools.gui.filter.*;
import com.insilico.application.insilicotools.gui.util.ImageUtil;

public class SubstructureFilterController extends GuiFilterController {
    public SubstructureFilterController(final String name, final TreeFilter treeFilter) {
        this(name, new SubstructureFilterState(), treeFilter);
    }

    public SubstructureFilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
        super(name, state, treeFilter);
    }

    protected Filter getFilter() {
        return new SubstructureFilter();
    }


    protected FilterGui getGuiComponent(FilterState state) {
        return new SubstructureFilterGui((SubstructureFilterState) state);
    }

    public FilterTransport getTransport() {
        return new SubstructureFilterTransport(state);
    }

    public String getToolTip() {
        String query = ((SubstructureFilterState) state).getQuery();
        if(query!=null) {
            return ImageUtil.smi2ImgUrl(query);
        }
        return "";
    }
}