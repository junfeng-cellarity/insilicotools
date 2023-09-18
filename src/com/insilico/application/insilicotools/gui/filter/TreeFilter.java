package com.insilico.application.insilicotools.gui.filter;

import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.MolPrinter;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.MolProperty;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.filter.alert.AlertFilterFactory;
import com.insilico.application.insilicotools.gui.filter.diverse.DiverseFilterFactory;
import com.insilico.application.insilicotools.gui.filter.dnd.DragAndDropLock;
import com.insilico.application.insilicotools.gui.filter.dnd.FilterPipelineTree;
import com.insilico.application.insilicotools.gui.filter.dnd.ScrollPaneTransferHandler;
import com.insilico.application.insilicotools.gui.filter.dnd.TableTransferHandler;
import com.insilico.application.insilicotools.gui.filter.protonate.ProtonatorFilterFactory;
import com.insilico.application.insilicotools.gui.table.MatrixCellTable;
import com.insilico.application.insilicotools.gui.table.MatrixMolTableModel;
import com.insilico.application.insilicotools.gui.filter.Neutralizer.NeutralizerFilterFactory;
import com.insilico.application.insilicotools.gui.filter.SDField.SDFieldFilterFactory;
import com.insilico.application.insilicotools.gui.filter.SaltStrip.SaltStripFilterFactory;
import com.insilico.application.insilicotools.gui.filter.deprotect.DeprotectFilterFactory;
import com.insilico.application.insilicotools.gui.filter.eganEgg.EganEggFilterFactory;
import com.insilico.application.insilicotools.gui.filter.fragment.FragmentFilterFactory;
import com.insilico.application.insilicotools.gui.filter.isotope.IsotopeFilterFactory;
import com.insilico.application.insilicotools.gui.filter.names.NameFilterFactory;
import com.insilico.application.insilicotools.gui.filter.property.PropertyFilterFactory;
import com.insilico.application.insilicotools.gui.filter.smarts.SmartsFilterFactory;
import com.insilico.application.insilicotools.gui.filter.substructure.SubstructureFilterFactory;
import com.insilico.application.insilicotools.gui.filter.uniq.UniqFilterFactory;
import com.insilico.application.insilicotools.gui.util.FileFunctor;
import com.insilico.application.insilicotools.gui.util.FileUtil;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolostream;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.tree.*;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.lang.String;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

public class TreeFilter extends JPanel {
	private final JTree tree;
	private final MatrixMolTableModel shoppingCartTableModel = new MatrixMolTableModel(new Vector<PropertyMolecule>());
	private final MatrixCellTable shoppingCartTable;
	private final MatrixMolTableModel displayTableModel = new MatrixMolTableModel(new Vector<PropertyMolecule>());
	private final MatrixCellTable displayTable = new MatrixCellTable(displayTableModel);
	private JPopupMenu popup = new JPopupMenu();
	private TreePath path;
	private final FilterFactory[] filters;
	private TitledBorder productsBorder;
	private Point editablePoint;
	protected EventListenerList shoppingCartListenerList = new EventListenerList();
	protected MoleculeSelectionEvent selectionEvent = null;
	private final List<FilterTreeNode> activeTopLevelFilterNodes = new ArrayList<FilterTreeNode>();
	private final List<FilterHolder> activeFilters = new ArrayList<FilterHolder>();
	private final HashMap<Object, FilterController> nodeToFilterControler = new HashMap<Object, FilterController>();
	public HashMap<FilterTreeNode, FilterHolder> filterNodeToHolder = new HashMap<FilterTreeNode, FilterHolder>();
	private JPanel shoppingCartPanel;
	private TitledBorder shoppingCartBorder;
    JProgressBar progressBar;
	JButton shoppingCartFiltersButton;

	public TreeFilter(PropertyMolecule[] initialMolecules) {
		super(new BorderLayout());
        progressBar = new JProgressBar(0,100);
        add(progressBar,BorderLayout.SOUTH);
		this.filters = new FilterFactory[]{
				new SubstructureFilterFactory(),
				new PropertyFilterFactory(),
				new EganEggFilterFactory(),
				new IsotopeFilterFactory(),
				new FragmentFilterFactory(),
				new SmartsFilterFactory(),
				new NameFilterFactory(),
				new UniqFilterFactory(),
				new SaltStripFilterFactory(),
				new NeutralizerFilterFactory(),
				new ProtonatorFilterFactory(),
				new DeprotectFilterFactory(),
				new SDFieldFilterFactory(initialMolecules),
                new AlertFilterFactory(),
				new DiverseFilterFactory()
        };

		tree = new FilterPipelineTree();

		ToolTipManager.sharedInstance().registerComponent(tree);

		shoppingCartTable = new MatrixCellTable(shoppingCartTableModel);

		ListSelectionListener selectionListener = new ListSelectionListener() {
			@Override
			public void valueChanged(ListSelectionEvent listSelectionEvent) {
				if (!listSelectionEvent.getValueIsAdjusting() && !tree.isSelectionEmpty()) {
					Object pathComponent = tree.getLastSelectedPathComponent();
					FilterTreeNode node = null;
					if (pathComponent instanceof FilterTreeNode) {
						node = (FilterTreeNode) pathComponent;
					} else if (pathComponent instanceof FilterResultTreeNode) {
						TreeNode previousNode = ((FilterResultTreeNode) pathComponent).getParent();
						if (previousNode instanceof FilterTreeNode) {
							node = (FilterTreeNode) previousNode;
						}
					}

					if (node == null) {
						return;
					}

				}
			}
		};
		shoppingCartTable.getSelectionModel().addListSelectionListener(selectionListener);
		shoppingCartTable.getColumnModel().getSelectionModel().addListSelectionListener(selectionListener);

		initialize(initialMolecules);

		if (initialMolecules != null) {
			displayTableModel.clear();
            System.out.println("loading molecules...");
			displayTableModel.addMolecules(Arrays.asList(initialMolecules));
            shoppingCartFiltersButton.setEnabled(displayTable.getSelectedMols().size()>0);
			updateBorder();
		}

	}

	public Dimension getPreferredSize() {
		return new Dimension(950, 660);
	}

	public void showPopup() {
		path = tree.getEditingPath();
		tree.cancelEditing();
		popup.show(tree, editablePoint.x, editablePoint.y);
	}

	public void addToShoppingCart() {
		PropertyMolecule[] molecules;

		Object o = tree.getSelectionPath().getLastPathComponent();
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

		shoppingCartTableModel.addMolecules(new ArrayList<PropertyMolecule>(Arrays.asList(molecules)));
		DragAndDropLock.setLocked(true);
		GhostGlassPane glassPane = (GhostGlassPane) (SwingUtilities.getRootPane(shoppingCartPanel).getGlassPane());
        try {
            BufferedImage itemPicture = getItemPicture(molecules);
            glassPane.setImage(itemPicture);
            glassPane.setVisible(true);
            glassPane.startAnimation(SwingUtilities.convertRectangle(shoppingCartTable, shoppingCartTable.getVisibleRect(), glassPane));
        } catch (MolFormatException e) {
            e.printStackTrace();
        }
    }

	private BufferedImage renderOffscreen(PropertyMolecule[] molecules) throws MolFormatException {
		if (molecules == null || molecules.length == 0) {
			return null;
		}
		int width = 150;
		int height = 150;
		BufferedImage image = new BufferedImage(width, height,
				BufferedImage.TYPE_INT_ARGB);
		Graphics2D imageGraphics = image.createGraphics();
		Rectangle r = new Rectangle(0, 0, width, height);

		Molecule molecule = MolImporter.importMol(ChemFunc.getMolString(molecules[0].getMol()));
		//MolImageSize imageSize = molecule.getImageSize(String.format("png:w%dh%d", width, height));
		MolPrinter molPrinter = new MolPrinter(molecule);
		molPrinter.setTransparent(true);
		molPrinter.setScale(molPrinter.maxScale(r));

		//molPrinter.setScale(imageSize.scale);
		molPrinter.paint(imageGraphics, r);
		for (int col = 0; col < width; col++) {
			for (int row = 0; row < height; row++) {
				if (image.getRGB(col, row) == Color.WHITE.getRGB()) image.setRGB(col, row, 0x0);
			}
		}
		return image;

	}

	public BufferedImage getItemPicture(PropertyMolecule[] molecules) throws MolFormatException {
		return renderOffscreen(molecules);
	}


	public void initialize(final PropertyMolecule[] initialMolecules) {
		if (initialMolecules != null) {
			reset(initialMolecules);
		} else {
			reset(new PropertyMolecule[0]);
		}

		ToolTipManager.sharedInstance().registerComponent(tree);

		//final TreeFilterRenderer renderer = new TreeFilterRenderer(tree, this);
		tree.setCellRenderer(new TreeFilterRenderer(tree, this));
		tree.setCellEditor(new TreeCellEditor() {
			private Object value;

			public Component getTreeCellEditorComponent(JTree tree, Object value, boolean isSelected, boolean expanded, boolean leaf, int row) {
				this.value = value;
				//return renderer.getTreeCellRendererComponent(tree, value, true, expanded, leaf, row, true);
				return new TreeFilterRenderer(tree, TreeFilter.this).getTreeCellRendererComponent(tree, value, true, expanded, leaf, row, true);
			}

			public void cancelCellEditing() {
			}

			public boolean stopCellEditing() {
				return true;
			}

			public Object getCellEditorValue() {
				return value;
			}

			public boolean isCellEditable(EventObject anEvent) {
				if (anEvent instanceof MouseEvent) {
					if (!((MouseEvent) anEvent).isPopupTrigger()) {
						editablePoint = ((MouseEvent) anEvent).getPoint();
						return true;
					}
				}
				return false;
			}

			public boolean shouldSelectCell(EventObject anEvent) {
				return true;
			}

			public void addCellEditorListener(CellEditorListener l) {
			}

			public void removeCellEditorListener(CellEditorListener l) {
			}
		});

		tree.setEditable(true);
		tree.setToggleClickCount(1);
		setPreferredSize(new Dimension(1000, 650));
		popup = new JPopupMenu();


		Arrays.sort(filters, new Comparator<Object>() {
			public int compare(Object o1, Object o2) {
				return o1.toString().compareToIgnoreCase(o2.toString());
			}
		});

		for (int i = 0; i < filters.length; i++) {
			final FilterFactory f = filters[i];
			JMenuItem item = new JMenuItem(f.getName());
			popup.add(item);
			item.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					runNewFilter(f.getInstance(TreeFilter.this));
				}
			});
		}

		final JPopupMenu bridgePopup = new JPopupMenu();

		JMenuItem renameNode = new JMenuItem("Rename");
		renameNode.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Object o = ((DefaultMutableTreeNode) path.getLastPathComponent()).getUserObject();
				String s = JOptionPane.showInputDialog("Rename", o == null ? "" : o.toString());
				if (s != null) {
					((DefaultMutableTreeNode) path.getLastPathComponent()).setUserObject(s);
				}
			}
		});
		bridgePopup.add(renameNode);

//        bridgePopup.add(deprotectRecursivelyNode);

		JMenuItem deleteNode = new JMenuItem("Delete");
		deleteNode.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				((DefaultTreeModel) tree.getModel()).removeNodeFromParent((MutableTreeNode) path.getLastPathComponent());
				FilterHolder holder = filterNodeToHolder.get(path.getLastPathComponent());
				activeFilters.remove(holder);  //hua need to do
				activeTopLevelFilterNodes.remove(path.getLastPathComponent());
			}
		});
		bridgePopup.add(deleteNode);


		tree.addMouseListener(new MouseAdapter() {
			public void mousePressed(MouseEvent e) {
				maybeShowPopup(e);
			}

			public void mouseReleased(MouseEvent e) {
				maybeShowPopup(e);
			}

			private void maybeShowPopup(MouseEvent e) {
				if (e.isPopupTrigger()) {
					JTree tree = (JTree) e.getComponent();
					path = tree.getPathForLocation(e.getX(), e.getY());
					if (path != null) {
						if (path.getLastPathComponent() instanceof FilterResultTreeNode) {
							tree.cancelEditing();
							popup.show(e.getComponent(), e.getX(), e.getY());
							tree.setSelectionPath(path);
						} else if (path.getLastPathComponent() != tree.getModel().getRoot()) {
							tree.cancelEditing();
							bridgePopup.show(e.getComponent(), e.getX(), e.getY());
							tree.setSelectionPath(path);
						}
					}
				}
			}
		});

		shoppingCartTable.setTransferHandler(new TableTransferHandler());
		JScrollPane tableScrollPane = new JScrollPane(shoppingCartTable);
		tableScrollPane.setTransferHandler(new ScrollPaneTransferHandler());

		shoppingCartPanel = new JPanel(new BorderLayout());
		shoppingCartPanel.add(tableScrollPane);
		shoppingCartBorder = new TitledBorder("Final Selection");
		shoppingCartPanel.setBorder(shoppingCartBorder);

		shoppingCartTableModel.addTableModelListener(new TableModelListener() {
			public void tableChanged(TableModelEvent e) {
				shoppingCartBorder.setTitle("Final Selection (" + shoppingCartTableModel.getPropertyMolecules().size() + ")");
				shoppingCartPanel.repaint();
				fireShoppingCartChanged(shoppingCartTableModel.getPropertyMolecules().toArray(new PropertyMolecule[shoppingCartTableModel.getPropertyMolecules().size()]), false);
			}
		});


		tree.getSelectionModel().addTreeSelectionListener(new TreeSelectionListener() {
			public void valueChanged(TreeSelectionEvent e) {
				if (!tree.isSelectionEmpty()) {
					Object o = tree.getLastSelectedPathComponent();
					if (o instanceof FilterResultTreeNode) {
						displayTableModel.clear();
						displayTableModel.addMolecules(Arrays.asList(((FilterResultTreeNode) o).getMolecules()));
						updateBorder();
						shoppingCartBorder.setTitle(String.format("Final Selection (%d)", shoppingCartTableModel.getPropertyMolecules().size()));
						shoppingCartPanel.repaint();
					}
				} else {
					displayTableModel.clear();
					updateBorder();
					shoppingCartBorder.setTitle(String.format("Final Selection (%d)", shoppingCartTableModel.getPropertyMolecules().size()));
					shoppingCartPanel.repaint();
				}
			}
		});


		JPanel treePanel = new JPanel(new BorderLayout());
		treePanel.add(new JScrollPane(tree));
		treePanel.setBorder(new TitledBorder("Filter Pipeline"));


		JSplitPane topSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treePanel, shoppingCartPanel);
		topSplitPane.setDividerLocation(400);

		JPanel displayTablePanel = new JPanel(new BorderLayout());
		productsBorder = new TitledBorder("Molecules in selected filter set");
		displayTablePanel.setBorder(productsBorder);
		displayTablePanel.add(new JScrollPane(displayTable));
		JToolBar displayTableToolbar = new JToolBar(JToolBar.HORIZONTAL);
		displayTableToolbar.setFloatable(false);
		displayTablePanel.add(displayTableToolbar, BorderLayout.SOUTH);

		shoppingCartFiltersButton = new JButton("Add selected to final selection", new ImageIcon(getClass().getClassLoader().getResource("shopping-cart.png")));
		shoppingCartFiltersButton.setToolTipText("Add selected to final selection");
//		shoppingCartFiltersButton.setBorder(null);
//		shoppingCartFiltersButton.setBorderPainted(false);
		shoppingCartFiltersButton.setEnabled(displayTable.getSelectedMols().size()>0);
		shoppingCartFiltersButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Vector<PropertyMolecule> selectedMols = displayTable.getSelectedMols();
				shoppingCartTableModel.addMolecules(selectedMols);
				DragAndDropLock.setLocked(true);
				GhostGlassPane glassPane = (GhostGlassPane) SwingUtilities.getRootPane(shoppingCartTable).getGlassPane();
				try {
					BufferedImage itemPicture = getItemPicture((PropertyMolecule[]) selectedMols.toArray(new PropertyMolecule[selectedMols.size()]));
					glassPane.setImage(itemPicture);
					glassPane.setVisible(true);
					glassPane.startAnimation(SwingUtilities.convertRectangle(shoppingCartTable, shoppingCartTable.getVisibleRect(), glassPane));
				} catch (MolFormatException e1) {
					e1.printStackTrace();
				}
			}
		});
		displayTableToolbar.add(shoppingCartFiltersButton);

		final JButton selectedAllButton = new JButton("Select all");
		selectedAllButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				displayTable.selectAll();
			}
		});
		displayTableToolbar.add(selectedAllButton);

		final JButton invertSelectionButton = new JButton("Invert Selection");
		invertSelectionButton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				displayTable.invertSelection();
			}
		});
		displayTableToolbar.add(invertSelectionButton);


		displayTable.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				if(shoppingCartFiltersButton!=null) {
					shoppingCartFiltersButton.setEnabled(displayTable.getNumSelected() > 0);
				}
			}
		});

		JSplitPane mainSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, topSplitPane, displayTablePanel);
		mainSplitPane.setDividerLocation(270);
		add(mainSplitPane, BorderLayout.CENTER);

		JToolBar filtertoolbar = new JToolBar("Filters", JToolBar.HORIZONTAL);
		filtertoolbar.setFloatable(false);

		filtertoolbar.addSeparator();

		JButton clearButton = new JButton("Clear", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Delete16.gif")));
		clearButton.setToolTipText("Clear");
		clearButton.setBorder(null);
		clearButton.setBorderPainted(false);
		clearButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				((DefaultTreeModel) tree.getModel()).setRoot(new FilterResultTreeNode(new FilterSection(((FilterResultTreeNode) tree.getModel().getRoot()).getMolecules(), "Molecule Set")));
				activeFilters.clear();
				activeTopLevelFilterNodes.clear();
			}
		});
		filtertoolbar.add(clearButton);
		treePanel.add(filtertoolbar, BorderLayout.SOUTH);




		JToolBar subReagentsToolbar = new JToolBar("Final Selection", JToolBar.HORIZONTAL);
		subReagentsToolbar.setFloatable(false);
		shoppingCartPanel.add(subReagentsToolbar, BorderLayout.SOUTH);

		final JButton sortButton = new JButton("Sort", new ImageIcon(Toolkit.getDefaultToolkit().createImage(getClass().getClassLoader().getResource("sorted.png"))));
		sortButton.setToolTipText("Sort");
		sortButton.setBorder(null);
		sortButton.setBorderPainted(false);
		sortButton.setEnabled(false);
		sortButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Component glassPane = getGlassPane();
				MouseAdapter glassPaneMouseListener = new MouseAdapter() {
				};
				glassPane.addMouseListener(glassPaneMouseListener);
				glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
				glassPane.setVisible(true);

				Collections.sort(shoppingCartTableModel.getPropertyMolecules(), new Comparator<PropertyMolecule>() {
					public int compare(PropertyMolecule o1, PropertyMolecule o2) {
						return new Double(o1.getMW()).compareTo(new Double(o2.getMW()));
					}
				});

				shoppingCartTableModel.fireTableDataChanged();

				glassPane.removeMouseListener(glassPaneMouseListener);
				glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
				glassPane.setVisible(false);
			}
		});

		subReagentsToolbar.add(sortButton);
		subReagentsToolbar.addSeparator();

		final JButton copyButton = new JButton("Copy", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Copy16.gif")));
		copyButton.setToolTipText("Copy");
		copyButton.setBorder(null);
		copyButton.setBorderPainted(false);
		copyButton.setEnabled(false);
		copyButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Component glassPane = getGlassPane();
				MouseAdapter glassPaneMouseListener = new MouseAdapter() {
				};
				glassPane.addMouseListener(glassPaneMouseListener);
				glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
				glassPane.setVisible(true);

				StringBuffer buffer = new StringBuffer();
				for (PropertyMolecule reagent : shoppingCartTableModel.getPropertyMolecules()) {
					if (reagent != null) {
						//buffer.append(reagent.getMolecule().toFormat("smiles"));
						buffer.append(reagent.getSmiles());
						buffer.append(" ");
						buffer.append(reagent.getName());
						buffer.append("\n");
					}
				}
				Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(buffer.toString()), null);

				glassPane.removeMouseListener(glassPaneMouseListener);
				glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
				glassPane.setVisible(false);
			}
		});

		subReagentsToolbar.add(copyButton);
		subReagentsToolbar.addSeparator();
		final JButton copyNamesButton = new JButton("Copy Names", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Copy16.gif")));
		copyNamesButton.setToolTipText("Copy Names");
		copyNamesButton.setBorder(null);
		copyNamesButton.setBorderPainted(false);
		copyNamesButton.setEnabled(false);
		copyNamesButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Component glassPane = getGlassPane();
				MouseAdapter glassPaneMouseListener = new MouseAdapter() {
				};
				glassPane.addMouseListener(glassPaneMouseListener);
				glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
				glassPane.setVisible(true);

				StringBuffer buffer = new StringBuffer();
				for (PropertyMolecule reagent : shoppingCartTableModel.getPropertyMolecules()) {
					if (reagent != null) {
						buffer.append(reagent.getName());
						buffer.append("\n");
					}
				}
				Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(buffer.toString()), null);

				glassPane.removeMouseListener(glassPaneMouseListener);
				glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
				glassPane.setVisible(false);
			}
		});

		subReagentsToolbar.add(copyNamesButton);
		subReagentsToolbar.addSeparator();

		final JButton clearTableButton = new JButton("Clear", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Delete16.gif")));
		clearTableButton.setToolTipText("Clear");
		clearTableButton.setBorder(null);
		clearTableButton.setBorderPainted(false);
		clearTableButton.setEnabled(false);
		clearTableButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				shoppingCartTableModel.clear();
				updateBorder();
			}
		});
		subReagentsToolbar.add(clearTableButton);
		subReagentsToolbar.addSeparator();

		final JButton saveSDFbutton = new JButton("Save sdf", new ImageIcon(getClass().getResource("/toolbarButtonGraphics/general/Save16.gif")));
		saveSDFbutton.setToolTipText("Save SDF");
		saveSDFbutton.setBorder(null);
		saveSDFbutton.setBorderPainted(false);
		saveSDFbutton.setEnabled(false);
		saveSDFbutton.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				final Component glassPane = getGlassPane();
				final MouseAdapter glassPaneMouseListener = new MouseAdapter() {
				};
				FileUtil.saveToFile(InSlilicoPanel.getInstance().getCurrentDirectory(), new FileNameExtensionFilter("SDF file","sdf"), new FileFunctor() {
					@Override
					public void execute(final File file) {
						SwingWorker sw = new SwingWorker() {
							@Override
							protected Object doInBackground() throws Exception {
								if(file!=null){
									File parent = file.getParentFile();
									if(parent!=null&&parent.exists()&&parent.isDirectory()){
										InSlilicoPanel.getInstance().setCurrentDirectory(parent);
									}
								}
								oemolostream ofs = new oemolostream();
								ofs.SetFormat(OEFormat.SDF);
								ofs.open(file.getAbsolutePath());
								int progress = 0;
								for(PropertyMolecule mol:getShoppingCart()){
									progress ++;
									Vector v = new Vector();
									v.add(String.format("Saving molecule No. %d",progress));
									v.add(100*progress/getShoppingCart().size());
									publish(v);
									oechem.OEWriteMolecule(ofs,mol.getMol());
								}
								ofs.close();
								return null;
							}

							@Override
							protected void process(List chunks) {
								glassPane.addMouseListener(glassPaneMouseListener);
								glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
								glassPane.setVisible(true);
							}

							@Override
							protected void done() {
								try {
									get();
									JOptionPane.showMessageDialog(TreeFilter.this, "SDF file saved.");
								} catch (Exception e1) {
									e1.printStackTrace();
									JOptionPane.showMessageDialog(TreeFilter.this,e1.getMessage());
								}finally {
									glassPane.removeMouseListener(glassPaneMouseListener);
									glassPane.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
									glassPane.setVisible(false);
								}
							}
						};
						sw.execute();

					}
				});
			}
		});
		subReagentsToolbar.add(saveSDFbutton);

		shoppingCartTableModel.addTableModelListener(new TableModelListener() {
			public void tableChanged(TableModelEvent e) {
				sortButton.setEnabled(shoppingCartTableModel.getRowCount() != 0);
				clearTableButton.setEnabled(shoppingCartTableModel.getRowCount() != 0);
				copyButton.setEnabled(shoppingCartTableModel.getRowCount() != 0);
				copyNamesButton.setEnabled(shoppingCartTableModel.getRowCount() != 0);
				saveSDFbutton.setEnabled(shoppingCartTableModel.getRowCount() != 0);
			}
		});
	}

	private void updateBorder() {
		productsBorder.setTitle("Molecules in selected filter set (" + displayTableModel.getPropertyMolecules().size() + ")");
		repaint();
	}

	public List<? extends PropertyMolecule> getShoppingCart() {
		return shoppingCartTableModel.getPropertyMolecules();
	}

    public void reset(PropertyMolecule[] molecules, String strName) {
        FilterResultTreeNode node = new FilterResultTreeNode(new FilterSection(molecules, strName));
        tree.cancelEditing();
        ((DefaultTreeModel) tree.getModel()).setRoot(node);
        tree.setSelectionInterval(0, 0);
        shoppingCartTableModel.clear();
        for (int i = 0; i < this.filters.length; i++) {
            FilterFactory filter = this.filters[i];
            if (filter instanceof SDFieldFilterFactory) {
                ((SDFieldFilterFactory) filter).updateProperties(molecules);
            }
        }
    }

	public void reset(PropertyMolecule[] molecules) {
        reset(molecules, "Molecule Set");
	}

    public void resetWithFilter(PropertyMolecule[] molecules, String strName) {
        FilterResultTreeNode node = new FilterResultTreeNode(new FilterSection(molecules, strName));
        tree.cancelEditing();
        ((DefaultTreeModel) tree.getModel()).setRoot(node);
        tree.setSelectionInterval(0, 0);
        shoppingCartTableModel.clear();
        for (int i = 0; i < this.filters.length; i++) {
            FilterFactory filter = this.filters[i];
            if (filter instanceof SDFieldFilterFactory) {
                ((SDFieldFilterFactory) filter).updateProperties(molecules);
            }
        }
    }

	private void runNewFilter(final FilterController filterController) {
		final FilterResultTreeNode parentNode = (FilterResultTreeNode) path.getLastPathComponent();
		final PropertyMolecule[] mols = parentNode.getMolecules();
        progressBar.setIndeterminate(true);
        final GhostGlassPane glassPane = (GhostGlassPane) (SwingUtilities.getRootPane(shoppingCartPanel).getGlassPane());
        glassPane.start();
        final SwingWorker<FilterResult,Object> sw = new SwingWorker<FilterResult, Object>() {
			@Override
			protected FilterResult doInBackground() throws Exception {
				return filterController.filter(new ProgressReporter() {
                    @Override
                    public void reportProgress(final String note, final int progress) {
                        SwingUtilities.invokeLater(new Runnable() {
                            @Override
                            public void run() {
                                progressBar.setString(note);
                                progressBar.setValue(progress);
                            }
                        });
                    }
                }, mols);
			}

			@Override
			protected void done() {
				try {
					FilterResult result = get();
					FilterTreeNode filterNode = new FilterTreeNode(filterController.getName(), filterController.getToolTip(), mols);
					if (result != null) {
						activeTopLevelFilterNodes.add(filterNode);
						nodeToFilterControler.put(filterNode, filterController);

						parentNode.add(filterNode);
						TreeNode[] nodes = parentNode.getPath();
						String strPath = "0";
						for (int i = 1; i < nodes.length; i++) {
							int indexChild = nodes[i - 1].getIndex(nodes[i]);
							strPath = strPath + Integer.toString(indexChild + 1);
						}
						FilterHolder holder = new FilterHolder(0, Integer.parseInt(strPath), filterController.getName(), filterController.getTransport());
						activeFilters.add(holder);
						filterNodeToHolder.put(filterNode, holder);
						FilterResultTreeNode firstNode = null;
						FilterSection[] sections = result.getSections();
						for (int i = 0; i < sections.length; i++) {
							FilterSection section = sections[i];
							FilterResultTreeNode sectionNode = new FilterResultTreeNode(section);
							if (firstNode == null) firstNode = sectionNode;
							filterNode.add(sectionNode);
						}

						((DefaultTreeModel) tree.getModel()).reload();
						if (firstNode != null) {
							TreePath firstPath = path.pathByAddingChild(filterNode).pathByAddingChild(firstNode);
							tree.scrollPathToVisible(firstPath);
							tree.getSelectionModel().setSelectionPath(firstPath);
						}
					}
				} catch (InterruptedException e) {
					e.printStackTrace();
				} catch (ExecutionException e) {
					e.printStackTrace();
				}finally {
                    progressBar.setIndeterminate(false);
					progressBar.setValue(0);
                    progressBar.setString("");
                    glassPane.stop();
				}
			}
		};
		sw.execute();
	}

	class FilterHolder {
		private final int id;
		private final String name;
		private final FilterTransport transport;
		private final int parentID;

		public FilterHolder(final int id, final String name, final FilterTransport transport) {
			this.id = id;
			this.name = name;
			this.transport = transport;
			this.parentID = 0;
		}

		public FilterHolder(final int id, final int pID, final String name, final FilterTransport transport) {
			this.id = id;
			this.name = name;
			this.transport = transport;
			this.parentID = pID;
		}

		public int getId() {
			return id;
		}

		public int getParentID() {
			return parentID;
		}

		public String getName() {
			return name;
		}

		public FilterTransport getTransport() {
			return transport;
		}

		public String toString() {
			return name;
		}
	}

	/*
	public void CreateSavedFilterTree() {

		try {
			if (this.filterSession > 0) {
				FilterHolder[] sessionFilters = getSavedSessionFilters();
				for (FilterHolder filterHolder : sessionFilters) {
					String name = filterHolder.getName();
					TreeModel treeModel = tree.getModel();
					TreeNode node = (TreeNode) treeModel.getRoot();
					runNewFilter(filterHolder.getTransport().getFilterController(name, TreeFilter.this), (FilterResultTreeNode) node, "");
				}
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}

	}
	*/


	public class FilterResultTreeNode extends DefaultMutableTreeNode {
		private PropertyMolecule[] molecules;

		public FilterResultTreeNode(FilterSection section) {
			super(section);
			this.molecules = section.getMolecules();
		}

        @Override
        public String getUserObject() {
            return ((FilterSection)super.getUserObject()).toString();
        }

        public PropertyMolecule[] getMolecules() {
			return molecules;
		}

		public void setMolecules(PropertyMolecule[] molecules) {
			this.molecules = molecules;
		}

		public String getName() {
			return getUserObject().toString();
		}

		public String toString() {
			return getUserObject() == null || molecules == null ? "" : getUserObject() + " (" + molecules.length + ")";
		}
	}

	public class FilterTreeNode extends DefaultMutableTreeNode {
		private final String tooltip;
		private final PropertyMolecule[] mols;

		public FilterTreeNode(String name, String tooltip, PropertyMolecule[] mols) {
			super(name);
			this.mols = mols;
			this.tooltip = tooltip;
		}

		public String getToolTip() {
			return tooltip;
		}

		public PropertyMolecule[] getMols() {
			return mols;
		}

	}

	public void addShoppingCartListener(MoleculeSelectionListener l) {
		shoppingCartListenerList.add(MoleculeSelectionListener.class, l);
	}

	public void removeShoppingCartListener(MoleculeSelectionListener l) {
		shoppingCartListenerList.remove(MoleculeSelectionListener.class, l);
	}

	protected void fireShoppingCartChanged(PropertyMolecule[] mols, boolean isValueAdjusting) {
		Object[] listeners = shoppingCartListenerList.getListenerList();
		for (int i = listeners.length - 2; i >= 0; i -= 2) {
			if (listeners[i] == MoleculeSelectionListener.class) {
				if (selectionEvent == null) selectionEvent = new MoleculeSelectionEvent(this, mols, isValueAdjusting);
				((MoleculeSelectionListener) listeners[i + 1]).moleculesSelected(selectionEvent);
			}
		}
		selectionEvent = null;
	}

	public Component getGlassPane() {
		Object o = getTopLevelAncestor();
		if (o instanceof JFrame) return ((JFrame) getTopLevelAncestor()).getGlassPane();
		return ((JDialog) getTopLevelAncestor()).getGlassPane();
	}
}