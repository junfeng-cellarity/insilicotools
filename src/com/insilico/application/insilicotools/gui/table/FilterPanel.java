package com.insilico.application.insilicotools.gui.table;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;
import com.jidesoft.swing.JideScrollPane;

import javax.swing.*;
import java.awt.*;
import java.util.Vector;

/**
 * Created by jfeng1 on 3/31/16.
 */
public class FilterPanel extends JPanel {
    TreeFilter filterTree;
    Vector<PropertyMolecule> molecules = new Vector<PropertyMolecule>();

    public FilterPanel() {
        super(new BorderLayout());
        filterTree = buildTree();
        add(new JideScrollPane(filterTree), BorderLayout.CENTER);
    }

    private TreeFilter buildTree(){
        return new TreeFilter(molecules.toArray(new PropertyMolecule[molecules.size()]));
    }

    public static void main(String[] args) {
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        FilterPanel p = new FilterPanel();
        f.getContentPane().add(p);
        f.setSize(new Dimension(1024,768));
        f.setVisible(true);

    }

}
