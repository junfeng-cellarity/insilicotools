package com.insilico.application.insilicotools.chart;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import org.jfree.data.xy.XYDataItem;

public class SelectionXYDataItem extends XYDataItem {

    private final Object key;
    private boolean selected;

    public SelectionXYDataItem(final double x, final double y, Object key) {
        super(x, y);
        this.key = key;
    }

    public Object getKey() {
        return key;
    }

    protected void setSelected(boolean selected) {
        this.selected = selected;
        if(key instanceof PropertyMolecule){
            ((PropertyMolecule)key).setIsSelected(selected);
        }
    }

    public boolean isSelected() {
        return selected;
    }

    public String toString() {
        return "x = " + super.getX() + ", y = " + super.getY() + ", key = " + key;
    }
}