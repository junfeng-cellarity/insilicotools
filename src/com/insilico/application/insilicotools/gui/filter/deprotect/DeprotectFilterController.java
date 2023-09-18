package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.insilico.application.insilicotools.gui.filter.*;


public class DeprotectFilterController extends GuiFilterController {
	public DeprotectFilterController(final String name, final TreeFilter treeFilter) {
		this(name, new DeprotectFilterState(), treeFilter);
	}

	protected DeprotectFilterController(String name, FilterState state, TreeFilter treeFilter) {
		super(name, state, treeFilter);
	}

	protected FilterGui getGuiComponent(FilterState state) {
		return new DeprotectFilterGui((DeprotectFilterState) state);  //To change body of implemented methods use File | Settings | File Templates.
	}

	protected Filter getFilter() {
		return new DeprotectFilter();  //To change body of implemented methods use File | Settings | File Templates.
	}

	protected boolean canUseDeprotectedSmiles() {
		return false;
	}

	public FilterTransport getTransport() {
		return new DeprotectFilterTransport(state);
	}
}
