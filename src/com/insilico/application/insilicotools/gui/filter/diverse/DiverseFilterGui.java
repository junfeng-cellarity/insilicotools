package com.insilico.application.insilicotools.gui.filter.diverse;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.widget.IntegerTextField;

import javax.swing.*;
import java.awt.*;

public class DiverseFilterGui extends FilterGui {
	private final IntegerTextField numOfCompoundsField = new IntegerTextField(5);

	private final DiverseFilterState state;
	public DiverseFilterGui(DiverseFilterState state1) {
		super(new FlowLayout(FlowLayout.CENTER));

		this.state = state1;

		add(new JLabel("Number of Diverse Compounds:"));
		add(numOfCompoundsField);
		if (state != null && state.isValidState() == null) {
			numOfCompoundsField.setText(""+state.getNumOfDiverseCompounds());
		}
	}

	public FilterState getState() {
		state.setNumOfDiverseCompounds(Integer.parseInt(numOfCompoundsField.getText()));
		return state;
	}
}
