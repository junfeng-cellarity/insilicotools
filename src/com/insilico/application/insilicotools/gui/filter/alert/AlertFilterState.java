package com.insilico.application.insilicotools.gui.filter.alert;

import com.insilico.application.insilicotools.gui.filter.FilterState;

import java.util.Vector;

public class AlertFilterState extends FilterState {
    private Vector<AlertRule> selectedRules;

    public String isValidState() {
        return selectedRules == null || selectedRules.size() == 0 ? "At least one alert rule must be selected" : null;
    }

    public Vector<AlertRule> getSelectedRules() {
        return selectedRules == null ? new Vector<AlertRule>() : selectedRules;
    }

    public void setSelectedRules(final Vector<AlertRule> selectedRules) {
        this.selectedRules = selectedRules;
    }
}