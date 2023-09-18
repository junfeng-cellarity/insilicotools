package com.insilico.application.insilicotools.gui.filter;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;

public abstract class FilterController {
	private String name;
	protected FilterState state;
	private Filter filter;
	protected final TreeFilter treeFilter;

	public FilterController(final String name, final FilterState state, final TreeFilter treeFilter) {
		this.name = name;
		this.state = state;
		this.treeFilter = treeFilter;
	}

	protected abstract Filter getFilter();

	public FilterTransport getTransport() {
		return null;
	}

	public FilterResult filter(final ProgressReporter progressReporter, final PropertyMolecule[] mols) throws Exception {

		if (filter == null) {
			filter = getFilter();
		}
		FilterResult result = filter.filter(progressReporter , mols, state);
		return result;
	}



	public String toString() {
		return getName();
	}

	public String getName() {
		return name;
	}

	public String getToolTip() {
		return null;
	}

	public void setName(String name) {
		this.name = name;
	}
}