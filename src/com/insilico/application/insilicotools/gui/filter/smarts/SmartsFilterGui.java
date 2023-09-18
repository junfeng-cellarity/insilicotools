package com.insilico.application.insilicotools.gui.filter.smarts;

import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class SmartsFilterGui extends FilterGui {
	private final JTextArea smartsTextArea;

	private final SmartsFilterState state;
	private final JCheckBox multipleMatchCB;

	public SmartsFilterGui(SmartsFilterState state1) {
		super(new BorderLayout());

		this.state = state1;

		smartsTextArea = new JTextArea(5, 40);
		smartsTextArea.setLineWrap(true);

		JPanel p = new JPanel(new BorderLayout());
		p.add(new JScrollPane(smartsTextArea));
		JPanel p2 = new JPanel();
		p2.add(new JLabel("Predefined Smarts:"));
		final JComboBox smartsCB = new JComboBox(SmartsPattern.getPredefinedSmartsPattern());
		smartsCB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent actionEvent) {
				SmartsPattern pattern = (SmartsPattern) smartsCB.getSelectedItem();
				smartsTextArea.setText(pattern.getSmarts());
			}
		});
		p2.add(smartsCB);

		multipleMatchCB = new JCheckBox("Multiple Match Only");
		multipleMatchCB.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent actionEvent) {
				state.setMoreThanOne(multipleMatchCB.isSelected());
			}
		});
		p2.add(multipleMatchCB);

		p.setBorder(new TitledBorder("Smarts"));

		add(p, BorderLayout.CENTER);
		add(p2, BorderLayout.SOUTH);
		if (state != null && state.isValidState() == null) {
			smartsTextArea.setText(state.getSmarts());
			multipleMatchCB.setSelected(state.isMoreThanOne());
		}
	}

	public FilterState getState() {
		state.setSmarts(smartsTextArea.getText());
		state.setMoreThanOne(multipleMatchCB.isSelected());
		return state;
	}
}
