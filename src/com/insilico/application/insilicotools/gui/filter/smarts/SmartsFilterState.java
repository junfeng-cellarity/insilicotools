package com.insilico.application.insilicotools.gui.filter.smarts;

import com.insilico.application.insilicotools.gui.filter.FilterState;

public class SmartsFilterState extends FilterState {
	private String smarts;
	private boolean moreThanOne = false;

	public String isValidState() {
		return smarts == null || smarts.trim().length() == 0 ? "You must enter a smarts string" : null;
	}

	public String getSmarts() {
		return smarts;
	}

	public void setSmarts(final String smarts) {
		this.smarts = smarts;
	}

	public boolean isMoreThanOne() {
		return moreThanOne;
	}

	public void setMoreThanOne(boolean moreThanOne) {
		this.moreThanOne = moreThanOne;
	}
}
