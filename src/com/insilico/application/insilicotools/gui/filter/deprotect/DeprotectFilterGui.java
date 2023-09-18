package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.google.common.collect.Lists;
import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;

import javax.swing.*;
import java.awt.*;


public class DeprotectFilterGui extends FilterGui {
	private final DeprotectFilterState state;
	java.util.List<JCheckBox> checkBoxes = Lists.newArrayList();
	DeprotectionGroup[] deprotectionGroups;


	public DeprotectFilterGui(DeprotectFilterState state) {
		super(new BorderLayout());
		try {
			deprotectionGroups = DeprotectionGroup.getDeprotectionGroups();
			JPanel checkBoxPanel = new JPanel();
			checkBoxPanel.setBorder(BorderFactory.createTitledBorder("Protecting Groups"));
			checkBoxPanel.setLayout(new BoxLayout(checkBoxPanel, BoxLayout.Y_AXIS));
			for (DeprotectionGroup group : deprotectionGroups) {
				JCheckBox checkBox = new JCheckBox(group.getName());
				checkBoxPanel.add(checkBox);
				checkBoxes.add(checkBox);
			}
			add(checkBoxPanel, BorderLayout.CENTER);

		} catch (Exception e) {
			e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
		}
		this.state = state;

	}

	private java.util.List<DeprotectionGroup> getGroups() {
		java.util.List<DeprotectionGroup> groupSmirks = Lists.newArrayList();
		if (deprotectionGroups != null) {
			for (JCheckBox groupCheckBox : checkBoxes) {
				if (groupCheckBox.isSelected()) {
					groupSmirks.add(this.deprotectionGroups[checkBoxes.indexOf(groupCheckBox)]);
				}
			}
		}
		return groupSmirks;
	}


	@Override
	public FilterState getState() {
		state.setDeprotectionGroups(getGroups());
		return state;  //To change body of implemented methods use File | Settings | File Templates.
	}
}
