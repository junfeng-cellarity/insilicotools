package com.insilico.application.insilicotools.gui.filter.names;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;

public class NameFilterGui extends FilterGui {
    private final JTextArea namesTextArea;

    private final NameFilterState state;

    public NameFilterGui(NameFilterState state) {
        super(new BorderLayout());

        this.state = state;

        namesTextArea = new JTextArea(15, 40);

        JPanel p = new JPanel(new BorderLayout());
        p.add(new JScrollPane(namesTextArea));
        p.setBorder(new TitledBorder("Names"));

        add(new JScrollPane(p), BorderLayout.CENTER);
    }

    public FilterState getState() {
        state.setNames(namesTextArea.getText());
        return state;
    }
}