package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.List;


public class DeprotectFilterState extends FilterState {
	List<DeprotectionGroup> deprotectionGroups;

	public DeprotectFilterState() {
	}

	public void setDeprotectionGroups(List<DeprotectionGroup> deprotectionGroups) {
		this.deprotectionGroups = deprotectionGroups;
	}

	public List<DeprotectionGroup> getDeprotectionGroups() {
		return deprotectionGroups;
	}

	@Override
	public String isValidState() {
		return (deprotectionGroups != null && deprotectionGroups.size() > 0) ? null : "At least one protection group need to be selected.";  //To change body of implemented methods use File | Settings | File Templates.
	}
}
