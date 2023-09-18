package com.insilico.application.insilicotools.gui.filter.dnd;

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.util.ArrayList;

import com.insilico.application.insilicotools.data.PropertyMolecule;

public class ArrayListTransferable implements Transferable {
    private final ArrayList<PropertyMolecule> data;
    private DataFlavor flavor;

    public ArrayListTransferable(ArrayList<PropertyMolecule> alist) {
        data = alist;

        try {
            flavor = new DataFlavor(DataFlavor.javaJVMLocalObjectMimeType + ";class=java.util.ArrayList");
        }
        catch(ClassNotFoundException err) {
            err.printStackTrace();
        }
    }

    public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException {
        if(!isDataFlavorSupported(flavor)) {
            throw new UnsupportedFlavorException(flavor);
        }
        return data;
    }

    public DataFlavor[] getTransferDataFlavors() {
        return new DataFlavor[]{flavor};
    }

    public boolean isDataFlavorSupported(DataFlavor flavor) {
        return flavor.equals(flavor);
    }
}