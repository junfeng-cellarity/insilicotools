package com.insilico.application.insilicotools.gui.filter;

import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;

import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;

public class TreeFilterRenderer extends DefaultTreeCellRenderer {
    private final Icon filterIcon = new ImageIcon(getClass().getClassLoader().getResource("filter.png"));
    private final Icon passIcon = new ImageIcon(getClass().getClassLoader().getResource("pass.png"));
    private final Icon failIcon = new ImageIcon(getClass().getClassLoader().getResource("fail.png"));
    private final Icon genericIcon = new ImageIcon(getClass().getClassLoader().getResource("generic.png"));
    private final Icon shoppingCartIcon = new ImageIcon(getClass().getClassLoader().getResource("shopping-cart.png"));
    private Icon rootIcon;

    private final DefaultTreeCellRenderer filterLabel = new DefaultTreeCellRenderer();
    private final Box filterPanel = new Box(BoxLayout.X_AXIS) {
        public String getToolTipText(MouseEvent event) {
            String tooltip = null;
            JComponent comp = filterPanel;
            Point p = event.getPoint();

            while (tooltip == null) {
                JComponent jchild = null;
                int count = comp.getComponentCount();
                for (int i = 0; i < count; i++) {
                    Component child = comp.getComponent(i);
                    if (child instanceof JComponent) {
                        Rectangle bounds = child.getBounds();
                        if (bounds.contains(p)) jchild = (JComponent) child;
                    }
                }
                if (jchild == null) break;
                comp = jchild;
                Point loc = comp.getLocation();
                p.translate(-loc.x, -loc.y);
                tooltip = comp.getToolTipText();
            }
            return tooltip;
        }
    };

    @SuppressWarnings("deprecation")
	public TreeFilterRenderer(JTree tree, final TreeFilter filter) {
        try {
            Molecule m = MolImporter.importMol("c1ccccc1");
            m.aromatize(false);
            rootIcon = new ImageIcon((Image) m.toObject("image:w18,h18"));
        }

        catch (Exception err) {
            err.printStackTrace();
        }

        filterPanel.add(filterLabel);
        JButton filterButton = new JButton(filterIcon);
        filterButton.setBorder(null);
        filterButton.setFocusPainted(false);
        filterButton.setBorderPainted(false);
        filterButton.setBackground(tree.getBackground());
        filterButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                filter.showPopup();
            }
        });
        filterButton.setToolTipText("Run a filter");
        filterPanel.add(filterButton);

        JButton shoppingCartButton = new JButton(shoppingCartIcon);
        shoppingCartButton.setBorder(null);
        shoppingCartButton.setFocusPainted(false);
        shoppingCartButton.setBorderPainted(false);
        shoppingCartButton.setBackground(tree.getBackground());
        shoppingCartButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                filter.addToShoppingCart();
            }
        });
        shoppingCartButton.setToolTipText("Move items to final selection");
        filterPanel.add(shoppingCartButton);
    }

    public Component getTreeCellRendererComponent(JTree tree, Object value, boolean sel, boolean expanded, boolean leaf, int row, boolean hasFocus) {
        setToolTipText(value instanceof TreeFilter.FilterTreeNode ? ((TreeFilter.FilterTreeNode) value).getToolTip() : null);
        if (!(value instanceof TreeFilter.FilterTreeNode)) {
            filterLabel.getTreeCellRendererComponent(tree, value, sel, expanded, leaf, row, hasFocus);

            Object userObject = ((DefaultMutableTreeNode) value).getUserObject();
            if (userObject instanceof FilterSection) {
                FilterSection section = (FilterSection) userObject;
                filterLabel.setIcon(section.getText().equals("Pass") || section.getState() == FilterSection.PASS ? passIcon : section.getText().equals("Fail") || section.getState() == FilterSection.FAIL ? failIcon : genericIcon);
                filterLabel.setToolTipText(section.getTooltip());
            } else {
                filterLabel.setIcon(rootIcon);
                filterLabel.setToolTipText(null);
            }

            filterPanel.revalidate();
            return filterPanel;
        } else {
            super.getTreeCellRendererComponent(tree, value, sel, expanded, leaf, row, hasFocus);
            setIcon(filterIcon);
        }
        return this;
    }
}