package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.filter.FilterTransport;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

public class DeprotectFilterTransport extends FilterTransport {
	protected DeprotectFilterTransport(FilterState state) {
		super(state);
	}

	@Override
	public FilterController getFilterController(String name, TreeFilter treeFilter) {
		return new DeprotectFilterController(name, state, treeFilter);
	}

}
