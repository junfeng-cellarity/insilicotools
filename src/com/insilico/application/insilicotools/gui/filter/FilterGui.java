package com.insilico.application.insilicotools.gui.filter;

import javax.swing.*;
import java.awt.*;

public abstract class FilterGui extends JPanel {
    public FilterGui() {
        super();
    }

    public FilterGui(LayoutManager layout) {
        super(layout);
    }

    public void onScreen() {

    }

    public void offScreen() {

    }

    public abstract FilterState getState();

    public String isValidState() {
        return getState().isValidState();
    }
}