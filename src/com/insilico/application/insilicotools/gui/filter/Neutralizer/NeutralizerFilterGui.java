package com.insilico.application.insilicotools.gui.filter.Neutralizer;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import java.awt.*;

public class NeutralizerFilterGui extends FilterGui {
    private final NeutralizerFilterState state;

    public NeutralizerFilterGui(NeutralizerFilterState state) {
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
