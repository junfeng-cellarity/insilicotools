package com.insilico.application.insilicotools.gui.filter.protonate;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.Neutralizer.NeutralizerFilterState;

import java.awt.*;

public class ProtonatorFilterGui extends FilterGui {
    private final ProtonatorFilterState state;

    public ProtonatorFilterGui(ProtonatorFilterState state) {
        super(new BorderLayout());
        this.state = state;
    }

    public FilterState getState() {
        return state;
    }

    public void addNotify() {
        super.addNotify();
    }
}
