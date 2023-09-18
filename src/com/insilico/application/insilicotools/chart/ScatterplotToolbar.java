package com.insilico.application.insilicotools.chart;

import javax.swing.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class ScatterplotToolbar extends JToolBar {

    private final Scatterplot plot;

    public ScatterplotToolbar(Scatterplot plot) {
        this.plot = plot;
    }

    public ScatterplotToolbar(int orientation, Scatterplot plot) {
        super(orientation);
        this.plot = plot;
        initialize();
    }

    public ScatterplotToolbar(String name, Scatterplot plot) {
        super(name);
        this.plot = plot;
        initialize();
    }

    public ScatterplotToolbar(String name, int orientation, Scatterplot plot) {
        super(name, orientation);
        this.plot = plot;
        initialize();
    }

    private void initialize() {
        JButton resetButton = new JButton("Reset");
        resetButton.setBorder(null);
		add(resetButton);
		resetButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                plot.reset();
            }
        });

        addSeparator();

        JButton selectAllButton = new JButton("Select All");
        selectAllButton.setBorder(null);
		add(selectAllButton);
		selectAllButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                plot.selectAll();
                plot.fireSelectionChanged();
                plot.repaint();
            }
        });

        addSeparator();

        JButton deSelectAllButton = new JButton("Deselect All");
        deSelectAllButton.setBorder(null);
		add(deSelectAllButton);
		deSelectAllButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                plot.deselectAll();
                plot.fireSelectionChanged();
                plot.repaint();
            }
        });

        addSeparator();

        JButton zoomInButton = new JButton("Zoom In");
        zoomInButton.setBorder(null);
		add(zoomInButton);
		zoomInButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                plot.zoomIn();
            }
        });

        addSeparator();

        JButton zoomOutButton = new JButton("Zoom Out");
        zoomOutButton.setBorder(null);
		add(zoomOutButton);
		zoomOutButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                plot.zoomOut();
            }
        });
    }
}