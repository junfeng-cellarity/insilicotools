package com.insilico.application.insilicotools.gui.util;

import org.jdesktop.swingx.JXTable;

import javax.swing.*;
import java.awt.*;

/**
 * Created by jfeng1 on 3/15/17.
 */
public class TableUtils {
    public static void scrollToVisible(JXTable table, int rowIndex) {
        table.getSelectionModel().setSelectionInterval(rowIndex, rowIndex);
        table.scrollRectToVisible(new Rectangle(table.getCellRect(rowIndex, 0, true)));
    }
}
