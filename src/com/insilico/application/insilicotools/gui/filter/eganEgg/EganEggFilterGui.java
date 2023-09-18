package com.insilico.application.insilicotools.gui.filter.eganEgg;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import javax.swing.*;
import java.awt.*;

public class EganEggFilterGui extends FilterGui {
    private final JRadioButton abs99RadioButton;
    private final JRadioButton abs95RadioButton;

    private final JRadioButton bbb99RadioButton;
    private final JRadioButton bbb90RadioButton;

    private final EganEggFilterState state;

    public EganEggFilterGui(EganEggFilterState state) {
        super(new BorderLayout());

        this.state = state;

        abs99RadioButton = new JRadioButton("99%", state.getAbsState() == EganEggFilterState.ABS_99);
        abs95RadioButton = new JRadioButton("95%", state.getAbsState() == EganEggFilterState.ABS_95);
        JRadioButton absIgnoreRadioButton = new JRadioButton("Ignore", state.getAbsState() == EganEggFilterState.IGNORE);

        bbb99RadioButton = new JRadioButton("99%", state.getBbbState() == EganEggFilterState.BBB_99);
        bbb90RadioButton = new JRadioButton("90%", state.getBbbState() == EganEggFilterState.BBB_90);
        JRadioButton bbbIgnoreRadioButton = new JRadioButton("Ignore", state.getBbbState() == EganEggFilterState.IGNORE);

        Box p = new Box(BoxLayout.Y_AXIS);

        JPanel absPanel = new JPanel();
        absPanel.add(new JLabel("Absorption Model"));
        absPanel.add(abs95RadioButton);
        absPanel.add(abs99RadioButton);
        absPanel.add(absIgnoreRadioButton);

        p.add(absPanel);

        JPanel bbbPanel = new JPanel();
        bbbPanel.add(new JLabel("Blood Brain Barrier Model"));
        bbbPanel.add(bbb90RadioButton);
        bbbPanel.add(bbb99RadioButton);
        bbbPanel.add(bbbIgnoreRadioButton);

        p.add(bbbPanel);

        ButtonGroup absGroup = new ButtonGroup();
        absGroup.add(abs99RadioButton);
        absGroup.add(abs95RadioButton);
        absGroup.add(absIgnoreRadioButton);

        ButtonGroup bbbGroup = new ButtonGroup();
        bbbGroup.add(bbb99RadioButton);
        bbbGroup.add(bbb90RadioButton);
        bbbGroup.add(bbbIgnoreRadioButton);

        add(p, BorderLayout.CENTER);
    }

    public FilterState getState() {
        state.setAbsState(abs99RadioButton.isSelected() ? EganEggFilterState.ABS_99 : abs95RadioButton.isSelected() ? EganEggFilterState.ABS_95 : EganEggFilterState.IGNORE);
        state.setBbbState(bbb99RadioButton.isSelected() ? EganEggFilterState.BBB_99 : bbb90RadioButton.isSelected() ? EganEggFilterState.BBB_90 : EganEggFilterState.IGNORE);
        return state;
    }
}