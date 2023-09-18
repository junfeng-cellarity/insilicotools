package com.insilico.application.insilicotools.chart;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class Series {

    private final String name;
    private final List points;

    public Series(String name) {
        this.name = name;
        this.points = new ArrayList();
    }

    public void addItem(SelectionXYDataItem item) {
        points.add(item);
    }

    public int size() {
        return points.size();
    }

    public boolean selectItem(Object key) {
        Iterator itr = points.iterator();
        while (itr.hasNext()) {
            SelectionXYDataItem item = (SelectionXYDataItem) itr.next();
            if (item.getKey().equals(key)) {
                item.setSelected(true);
                return true;
            }
        }
        return false;
    }

    public SelectionXYDataItem getItem(int item) {
        return (SelectionXYDataItem) points.get(item);
    }

    public void clearSelection() {
        Iterator itr = points.iterator();
        while (itr.hasNext()) {
            ((SelectionXYDataItem) itr.next()).setSelected(false);
        }
    }

    public void selectAll() {
        Iterator itr = points.iterator();
        while (itr.hasNext()) {
            ((SelectionXYDataItem) itr.next()).setSelected(true);
        }
    }

    public void deselectAll() {
        Iterator itr = points.iterator();
        while (itr.hasNext()) {
            ((SelectionXYDataItem) itr.next()).setSelected(false);
        }
    }

    public SelectionXYDataItem[] getSelectedPoints() {
        List selectedPoints = new ArrayList();
        Iterator itr = points.iterator();
        while (itr.hasNext()) {
            SelectionXYDataItem item = (SelectionXYDataItem) itr.next();
            if (item.isSelected()) {
                selectedPoints.add(item);
            }
        }
        return (SelectionXYDataItem[]) selectedPoints.toArray(new SelectionXYDataItem[0]);
    }

    public String getName() {
        return name;
    }
}