package com.insilico.application.insilicotools.gui.filter.deprotect;

import com.insilico.application.insilicotools.gui.filter.FilterController;
import com.insilico.application.insilicotools.gui.filter.FilterFactory;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;


public class DeprotectFilterFactory extends FilterFactory {
	public DeprotectFilterFactory() {
		super("Deprotect");
	}

	@Override
	public FilterController getInstance(TreeFilter treeFilter) {
		return new DeprotectFilterController(getName(), treeFilter);  //To change body of implemented methods use File | Settings | File Templates.
	}
}
