package com.insilico.application.insilicotools.gui.filter.fragment;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: liuh
 * Date: 6/26/12
 * Time: 4:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class FragmentFilterGui extends FilterGui {
    private FragmentFilterState state;

    public FragmentFilterGui(FragmentFilterState state) {
        super(new BorderLayout());
        this.state = state;
    }
    public FilterState getState() {
        return state;
    }
}
