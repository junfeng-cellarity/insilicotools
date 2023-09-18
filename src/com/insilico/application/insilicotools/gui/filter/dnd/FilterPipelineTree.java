package com.insilico.application.insilicotools.gui.filter.dnd;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.filter.GhostGlassPane;
import com.insilico.application.insilicotools.gui.filter.TreeFilter;

import javax.swing.*;
import java.awt.*;
import java.awt.dnd.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;

public class FilterPipelineTree extends JTree implements DragGestureListener, DragSourceMotionListener, DragSourceListener {
	private static final int ITEM_WIDTH = 300;
	private static final int ITEM_HEIGHT = 300;
	private Point origin;

	public FilterPipelineTree() {
		DragSource dragSource = new DragSource();
		dragSource.createDefaultDragGestureRecognizer(this, DnDConstants.ACTION_COPY, this);
		dragSource.addDragSourceMotionListener(this);
	}

	public static int getItemWidth() {
		return ITEM_WIDTH;
	}

	public static int getItemHeight() {
		return ITEM_HEIGHT;
	}

//	private BufferedImage renderOffscreen(PropertyMolecule[] molecules) {
//		Vector<PropertyMolecule> molsForImage = new Vector<PropertyMolecule>();
//		for (int i = 0; i < Math.min(molecules.length, 50); i++){
//			molsForImage.add(molecules[i]);
//		}
//		return ChemFunc.getMolGridImage(molsForImage,5,5);
//	}

//	public BufferedImage getItemPicture(PropertyMolecule[] molecules) {
//		return renderOffscreen(molecules);
//	}

	public void dragGestureRecognized(DragGestureEvent dge) {
		if (DragAndDropLock.isLocked()) {
			DragAndDropLock.setDragAndDropStarted(false);
			return;
		}
		DragAndDropLock.setLocked(true);
		DragAndDropLock.setDragAndDropStarted(true);

		origin = (Point) dge.getDragOrigin().clone();

		PropertyMolecule[] molecules;

		Object o = getSelectionPath().getLastPathComponent();
		if (o instanceof TreeFilter.FilterResultTreeNode) {
			molecules = ((TreeFilter.FilterResultTreeNode) o).getMolecules();
		} else {
			TreeFilter.FilterTreeNode node = (TreeFilter.FilterTreeNode) o;
			ArrayList<PropertyMolecule> mols = new ArrayList<PropertyMolecule>();
			for (int i = 0; i < node.getChildCount(); i++) {
				mols.addAll(Arrays.asList(((TreeFilter.FilterResultTreeNode) node.getChildAt(i)).getMolecules()));
			}
			molecules = mols.toArray(new PropertyMolecule[mols.size()]);
		}

		dge.startDrag(Cursor.getDefaultCursor(), new ArrayListTransferable(new ArrayList<PropertyMolecule>(Arrays.asList(molecules))), this);



		GhostGlassPane glassPane = (GhostGlassPane) SwingUtilities.getRootPane(this).getGlassPane();
		glassPane.setVisible(true);

		Point p = (Point) dge.getDragOrigin().clone();
		SwingUtilities.convertPointToScreen(p, this);
		SwingUtilities.convertPointFromScreen(p, glassPane);

		glassPane.setPoint(p);
//		glassPane.setImage(getItemPicture(molecules), ITEM_WIDTH);
		glassPane.repaint();
	}

	public void dragMouseMoved(DragSourceDragEvent dsde) {
		if (!DragAndDropLock.isDragAndDropStarted()) {
			return;
		}

		GhostGlassPane glassPane = (GhostGlassPane) SwingUtilities.getRootPane(this).getGlassPane();

		Point p = (Point) dsde.getLocation().clone();
		SwingUtilities.convertPointFromScreen(p, glassPane);
		glassPane.setPoint(p);

		glassPane.repaint(glassPane.getRepaintRect());
	}

	public void dragEnter(DragSourceDragEvent dsde) {
	}

	public void dragOver(DragSourceDragEvent dsde) {
	}

	public void dropActionChanged(DragSourceDragEvent dsde) {
	}

	public void dragExit(DragSourceEvent dse) {
	}

	public void dragDropEnd(DragSourceDropEvent dsde) {
		if (!DragAndDropLock.isDragAndDropStarted()) {
			return;
		}
		DragAndDropLock.setDragAndDropStarted(false);

		GhostGlassPane glassPane = (GhostGlassPane) SwingUtilities.getRootPane(this).getGlassPane();

		Point p = (Point) dsde.getLocation().clone();
		SwingUtilities.convertPointFromScreen(p, glassPane);

		if (!dsde.getDropSuccess()) {
			SwingUtilities.convertPointToScreen(origin, this);
			SwingUtilities.convertPointFromScreen(origin, glassPane);

			Timer backTimer = new Timer(1000 / 60, new TravelBackToOrigin(glassPane, p, origin));
			backTimer.start();
		} else {
			glassPane.setPoint(p);
			glassPane.repaint(glassPane.getRepaintRect());
		}
	}

	private class TravelBackToOrigin implements ActionListener {
		private boolean isInitialized;
		private long start;

		private Point startPoint;
		private Point endPoint;
		private GhostGlassPane glassPane;

		private static final double INITIAL_SPEED = 500.0;
		private static final double INITIAL_ACCELERATION = 6000.0;

		private TravelBackToOrigin(GhostGlassPane glassPane, Point start, Point end) {
			this.glassPane = glassPane;
			this.startPoint = start;
			this.endPoint = end;
			isInitialized = false;
		}

		public void actionPerformed(ActionEvent e) {
			if (!isInitialized) {
				isInitialized = true;
				start = System.currentTimeMillis();
			}

			long elapsed = System.currentTimeMillis() - start;
			double time = elapsed / 1000.0;

			double a = (endPoint.y - startPoint.y) / (double) (endPoint.x - startPoint.x);
			double b = endPoint.y - a * endPoint.x;

			int travelX = (int) (INITIAL_ACCELERATION * time * time / 2.0 + INITIAL_SPEED * time);
			if (startPoint.x > endPoint.x) {
				travelX = -travelX;
			}

			int travelY = (int) ((startPoint.x + travelX) * a + b);
			int distanceX = Math.abs(startPoint.x - endPoint.x);

			if (Math.abs(travelX) >= distanceX) {
				((Timer) e.getSource()).stop();

				glassPane.setPoint(endPoint);
				glassPane.repaint(glassPane.getRepaintRect());

				SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						glassPane.setImage(null);
						glassPane.setVisible(false);
					}
				});
				DragAndDropLock.setLocked(false);

				return;
			}

			glassPane.setPoint(new Point(startPoint.x + travelX, travelY));

			glassPane.repaint(glassPane.getRepaintRect());
		}
	}
}