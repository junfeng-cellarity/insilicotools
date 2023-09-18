package com.insilico.application.insilicotools.gui.filter.uniq;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.awt.*;

public class UniqFilterGui extends FilterGui {
    private final UniqFilterState state;

    public UniqFilterGui(UniqFilterState state) {
        super(new BorderLayout());
        this.state = state;
    }

    public FilterState getState() {
        return state;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void addNotify() {
        super.addNotify();
    }

}
