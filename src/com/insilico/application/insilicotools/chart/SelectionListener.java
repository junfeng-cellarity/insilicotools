package com.insilico.application.insilicotools.chart;

import java.util.EventListener;

public interface SelectionListener extends EventListener {

    public void selectionChanged(SelectionEvent e);
}