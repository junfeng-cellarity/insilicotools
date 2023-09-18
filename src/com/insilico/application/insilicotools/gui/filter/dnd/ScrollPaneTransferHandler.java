package com.insilico.application.insilicotools.gui.filter.dnd;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.filter.GhostGlassPane;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;

import javax.swing.*;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.IOException;
import java.util.ArrayList;

public class ScrollPaneTransferHandler extends TransferHandler {
    private DataFlavor flavor;

    public ScrollPaneTransferHandler() {
        try {
            flavor = new DataFlavor(DataFlavor.javaJVMLocalObjectMimeType + ";class=java.util.ArrayList");
        }
        catch(ClassNotFoundException err) {
            err.printStackTrace();
        }
    }

    public int getSourceActions(JComponent c) {
        return NONE;
    }

    public boolean importData(JComponent c, Transferable t) {
        if(canImport(c, t.getTransferDataFlavors())) {
            try {
                @SuppressWarnings("unchecked")
				ArrayList<PropertyMolecule> data = (ArrayList<PropertyMolecule>) t.getTransferData(flavor);
                if(data != null) {
                    MatrixMolTableModel tableModel = (MatrixMolTableModel) ((JTable) ((JScrollPane) c).getViewport().getView()).getModel();
                    tableModel.addMolecules(data);
                }

                GhostGlassPane glassPane = (GhostGlassPane) SwingUtilities.getRootPane(c).getGlassPane();
                glassPane.startAnimation(SwingUtilities.convertRectangle(c, c.getVisibleRect(), glassPane));
                return true;
            }
            catch(UnsupportedFlavorException ufe) {
            }
            catch(IOException ioe) {
            }
        }

        return false;
    }

    public boolean canImport(JComponent c, DataFlavor[] flavors) {
        for(int i = 0; i < flavors.length; i++) {
            if(flavor.equals(flavors[i])) {
                return true;
            }
        }
        return false;
    }
}