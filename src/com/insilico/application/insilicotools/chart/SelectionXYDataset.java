package com.insilico.application.insilicotools.chart;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jfree.data.xy.AbstractXYDataset;

public class SelectionXYDataset extends AbstractXYDataset {

    private final List data;

    public SelectionXYDataset() {
        this.data = new ArrayList();
    }

    public int addSeries(String name) {
        data.add(new Series(name));
        fireDatasetChanged();
        return data.size() - 1;
    }

    public void deleteSeries(int series) {
        data.remove(series);
        fireDatasetChanged();
    }

    public boolean selectPoint(int series, Object key) {
        return ((Series) data.get(series)).selectItem(key);
    }

    public void addItem(int series, SelectionXYDataItem item) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        ((Series) data.get(series)).addItem(item);
        fireDatasetChanged();
    }

    public int getItemCount(int series) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        return ((Series) data.get(series)).size();
    }

    public Number getX(int series, int item) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        return ((Series) data.get(series)).getItem(item).getX();
    }

    public Number getY(int series, int item) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        return ((Series) data.get(series)).getItem(item).getY();
    }

    public void clearSelection() {
        for (int i = 0; i < data.size(); i++) {
            ((Series) data.get(i)).clearSelection();
        }
    }

    public void selectAll() {
        for (int i = 0; i < data.size(); i++) {
            selectAll(i);
        }
    }

    public void selectAll(int series) {
        ((Series) data.get(series)).selectAll();
    }

    public SelectionXYDataItem getItem(int series, int item) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        return ((Series) data.get(series)).getItem(item);
    }

    public int getSeriesCount() {
        return data.size();
    }

	public Comparable getSeriesKey(int series) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        return ((Series) data.get(series)).getName();
	}

    
    public String getSeriesName(int series) {
        if (series >= getSeriesCount()) {
            throw new ScatterplotRuntimeException("Undefined series");
        }
        return ((Series) data.get(series)).getName();
    }

    public SelectionXYDataItem[] getSelectedPoints() {
        List selectedPoints = new ArrayList();
        for (int i = 0; i < data.size(); i++) {
            selectedPoints.addAll(Arrays.asList(((Series) data.get(i)).getSelectedPoints()));
        }
        return (SelectionXYDataItem[]) selectedPoints.toArray(new SelectionXYDataItem[0]);
    }

    public void clear() {
        data.clear();
        fireDatasetChanged();
    }

    public void deselectAll(int series) {
        ((Series) data.get(series)).deselectAll();
    }

    public void deselectAll() {
        for (int i = 0; i < data.size(); i++) {
            deselectAll(i);
        }
    }

}