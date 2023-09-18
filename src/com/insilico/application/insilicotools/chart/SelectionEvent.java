package com.insilico.application.insilicotools.chart;

import java.util.EventObject;

public class SelectionEvent extends EventObject {

    private final Object[] keys;

    public SelectionEvent(Object source, SelectionXYDataItem[] selectedPoints) {
        super(source);
        this.keys = new Object[selectedPoints.length];
        for (int i = 0; i < selectedPoints.length; i++) {
            keys[i] = selectedPoints[i].getKey();
        }
    }

    public Object[] getSelectedKeys() {
        return keys;
    }
}