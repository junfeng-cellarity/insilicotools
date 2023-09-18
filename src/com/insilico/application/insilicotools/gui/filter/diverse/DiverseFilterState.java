package com.insilico.application.insilicotools.gui.filter.diverse;

import com.insilico.application.insilicotools.gui.filter.FilterState;

public class DiverseFilterState extends FilterState {
	private int numOfDiverseCompounds;

	public String isValidState() {
		if(numOfDiverseCompounds<=0){
			return "Number of Compounds required.";
		}else{
			return null;
		}
	}

	public int getNumOfDiverseCompounds() {
		return numOfDiverseCompounds;
	}

	public void setNumOfDiverseCompounds(int numOfDiverseCompounds) {
		this.numOfDiverseCompounds = numOfDiverseCompounds;
	}
}
