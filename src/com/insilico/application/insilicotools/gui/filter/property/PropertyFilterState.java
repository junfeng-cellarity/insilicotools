package com.insilico.application.insilicotools.gui.filter.property;

import com.insilico.application.insilicotools.gui.filter.FilterState;

public class PropertyFilterState extends FilterState {
    private PropertyCondition[] conditions;

    public PropertyCondition[] getConditions() {
        return conditions;
    }

    public void setConditions(final PropertyCondition[] conditions) {
        this.conditions = conditions;
    }

    public String isValidState() {
        return conditions.length == 0 ? "At least one condition must be specified" : null;
    }
}