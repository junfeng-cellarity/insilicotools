package com.insilico.application.insilicotools.gui.widget;

import javax.swing.*;
import javax.swing.table.DefaultTableCellRenderer;
import java.awt.*;

public class StarRaterCellRender extends DefaultTableCellRenderer{
    public StarRaterCellRender() {
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        if(value!=null){
            int rank = (Integer)value;
            return new StarRater(5,(float)rank);
        }else{
            return new StarRater(5,0);
        }
    }

}
