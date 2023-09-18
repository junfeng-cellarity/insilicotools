package com.insilico.application.insilicotools.gui.filter.substructure;

import chemaxon.formats.MolExporter;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.gui.filter.FilterGui;
import com.insilico.application.insilicotools.gui.filter.FilterState;
import com.insilico.application.insilicotools.gui.util.MarvinFactory;

import java.awt.*;
import java.awt.event.KeyEvent;
import java.io.IOException;
import java.util.HashMap;

public class SubstructureFilterGui extends FilterGui {
    private final MSketchPane sketcher = MarvinFactory.getCompoundSketcher();
    private final SubstructureFilterState state;
    protected KeyEventPostProcessor keyEventProcessor;

    public SubstructureFilterGui(SubstructureFilterState state) {
        this.state = state;
        add(sketcher, BorderLayout.CENTER);
        if(state.getQuery()!=null) {
            sketcher.setMol(state.getQuery());
        }

        keyEventProcessor = new KeyEventPostProcessor() {
            private HashMap<Object, KeyEvent> unconsumedEvents = new HashMap<Object, KeyEvent>();

            public boolean postProcessKeyEvent(KeyEvent e) {
                if (e.getModifiers() > 0) return false;

                maybeProcess(e);
                return false;
            }

            private void maybeProcess(KeyEvent e) {
                Object source = e.getSource();
                KeyEvent lastEvt = unconsumedEvents.remove(source);
                switch (e.getID()) {
                    case KeyEvent.KEY_PRESSED:
                        if (!e.isConsumed()) unconsumedEvents.put(source, e);
                        break;
                    case KeyEvent.KEY_TYPED:
                        if (!e.isConsumed() && lastEvt != null) {
//                            sketcher.keyPressed(lastEvt);
//                            sketcher.keyTyped(e);
                            unconsumedEvents.put(source, e);
                        }
                        break;
                    case KeyEvent.KEY_RELEASED:
//                        if (!e.isConsumed() && lastEvt != null) sketcher.keyReleased(e);
                        break;
                }
            }
        };
    }

    public void onScreen() {
        KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventPostProcessor(keyEventProcessor);
    }

    public void offScreen() {
        KeyboardFocusManager.getCurrentKeyboardFocusManager().removeKeyEventPostProcessor(keyEventProcessor);
    }

    public FilterState getState() {
        try {
            Molecule mol = sketcher.getMol();
            String smarts = MolExporter.exportToFormat(mol, "smarts:ah");
            state.setQuery(smarts);
            return state;
        } catch (IOException e) {
            state.setQuery(null);
            e.printStackTrace();
        }
        return state;
    }
}