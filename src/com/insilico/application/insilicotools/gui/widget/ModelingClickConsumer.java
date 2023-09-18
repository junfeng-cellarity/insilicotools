package com.insilico.application.insilicotools.gui.widget;

import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.CompoundInputDialog;
import com.insilico.application.insilicotools.gui.DesignProgressMonitor;
import com.insilico.application.insilicotools.gui.dialog.MolBuildingDialog;
import com.insilico.application.insilicotools.gui.modeling.ModelingException;
import com.insilico.application.insilicotools.gui.util.ImageUtil;
import com.insilico.application.insilicotools.gui.util.JyMolUtilities;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.insilico.application.insilicotools.util.Protonator;
import com.schrodinger.jymol.JyMol;
import com.schrodinger.jymol.JyMolClickConsumer;
import openeye.oechem.*;
import org.apache.xmlrpc.XmlRpcException;

import javax.swing.*;
import javax.swing.event.EventListenerList;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.net.MalformedURLException;
import java.util.*;
import java.util.List;
import java.util.concurrent.ExecutionException;

public class ModelingClickConsumer extends JyMolClickConsumer implements JymolClickedEventListener{
	final protected JyMol jymol;
	MolViewer3D molViewer3D;
	Hashtable<String, String> psName2dName = new Hashtable();
	int numDistances = 0;
	int numAngles = 0;
	int numDihedrals = 0;
	//
	//Define a torsion for rotation: A(torsionStart)-----
	//                                                   \
	//                                                    \
	//                                                     B(bondStart)----C(bondEnd)
	//                                                                     \
	//                                                                      \
	//                                                                       D(torsionEnd)
	JymolSelection torsionStart;
	JymolSelection bondStart;
	JymolSelection bondEnd;
	JymolSelection torsionEnd;

	MolObject molObject;
	MolObject fragObject;

	Vector<MolObject> molObjects;
	HashMap<String,MolObject> molObjectHash;
	MolObject receptorObj;
    boolean showHBond = true;
    MolObject currentHbondMolObj = null;
    String receptorName;
    CompoundInputDialog superimposeDialog;


	Vector<JymolSelection> selectionBuffer = new Vector<JymolSelection>();
	protected JPopupMenu menu = new JPopupMenu();
	public final static int DISTANCE_MODE = 1;
	public final static int ANGLE_MODE = 2;
	public final static int DIHEDRAL_MODE = 3;
	public final static int BUILDING_MODE = 5;
	public final static int SELECTION_MODE = 0;
	public final static int TORSION_SELECTION_MODE = 4;
    public final static int LABELING_MODE = 14;

	public final static int BUILDING_SUBMODE_VIEW = 6;
	public final static int BUILDING_SUBMODE_ATTACH = 7;
	public final static int BUILDING_SUBMODE_MODIFY = 8;
	public final static int BUILDING_SUBMODE_LINK = 9;
	public final static int BUILDING_SUBMODE_DELETE = 10;
	public final static int BUILDING_SUBMODE_MODIFY_BOND = 11;
	public final static int BUILDING_SUBMODE_DELETE_BOND = 12;
	public final static int BUILDING_SUBMODE_ROTATE_BOND = 13;

	protected int currentMode = SELECTION_MODE;
	int building_mode = BUILDING_SUBMODE_VIEW;

	MorphingFragment morphingFragment;
	MorphingElement morphingElement;
	int morphing_linker_length = -1;
	String morphing_bond_type = "single";
	String receptorStyle = JyMolUtilities.LINE_STYLE;
	String ligandStyle = JyMolUtilities.STICK_STYLE;
	String receptorColor = "green";

    boolean zoom = true;
	boolean receptorIsShowing = false;

	int hydrogen_option = SHOW_NO_HYDROGEN;
	int receptor_hydrogen_option = SHOW_NO_HYDROGEN;
	public final static int SHOW_NO_HYDROGEN = 0;
	public final static int SHOW_POLAR_HYDROGEN = 1;
	public final static int SHOW_ALL_HYDROGEN = 2;

	JRadioButtonMenuItem buildMenuItem;

	JMenu selectedMoleculeMenu;
    
    Vector<String> labelBuffer = new Vector<String>();

	Stack<Integer> matchList = new Stack<Integer>();

	MolBuildingDialog molBuildingDialog;

	Stack<MolObject> undoStack = new Stack<MolObject>();
	Stack<MolObject> redoStack = new Stack<MolObject>();

	public void undo(){
		if(!undoStack.isEmpty()){
			MolObject molObj = undoStack.pop();
			MolObject oldObj = molObjectHash.get(molObj.getName());
			if(oldObj!=null) {
				updateMolObject(molObj, false);
				jymol.firePropertyChange(molObj.getName(),true,false);
				redoStack.push(oldObj);
			}
		}
	}

	public void redo(){
		if(!redoStack.isEmpty()){
			MolObject molObj = redoStack.pop();
			MolObject oldObj = molObjectHash.get(molObj.getName());
			if(oldObj!=null) {
				updateMolObject(molObj, false);
                jymol.firePropertyChange(molObj.getName(),true,false);
				undoStack.push(oldObj);
			}
		}
	}

	protected EventListenerList listenerList = new EventListenerList();

	public void addJymolClickedEventListener(JymolClickedEventListener listener) {
		listenerList.add(JymolClickedEventListener.class, listener);
	}

	public void removeJymolClickedEventListener(JymolClickedEventListener listener) {
		listenerList.remove(JymolClickedEventListener.class, listener);
	}

	void fireJymolClickedEvent(JymolClickedEvent event) {
		Object[] listeners = listenerList.getListenerList();
		for (int i = 0; i < listeners.length; i += 2) {
			if (listeners[i] == JymolClickedEventListener.class) {
				((JymolClickedEventListener) listeners[i + 1]).jymolClicked(event);
			}
		}
	}

    public MolViewer3D getMolViewer3D() {
        return molViewer3D;
    }

    public Vector getMatchList() {
		return matchList;
	}

	public void addSelection(String molName, Vector<Integer> selection) {
        clearSelection();
		for (int i = 0; i < matchList.size(); i++) {
			Integer sel = (Integer) matchList.elementAt(i);
			JyMolUtilities.label(jymol, sel.toString(), "");
			jymol.cmd.delete(sel.toString());
		}
		matchList.clear();

		matchList.addAll(selection);
		for (int i = 0; i < matchList.size(); i++) {
			Integer id = matchList.get(i);
			String selectionName = id.toString();
			String index = new Integer(i + 1).toString();
			int[] selected = new int[1];
			selected[0] = id.intValue();
			jymol.cmd.selectList(selectionName, molName, selected);
			JyMolUtilities.label(jymol, selectionName, index);
		}
		if(molObjectHash.containsKey(molName)){
			molObject = molObjectHash.get(molName);
		}
		fireJymolClickedEvent(new JymolClickedEvent(this));
	}

	public void clearSelection() {
		molObject = null;
		for (int i = 0; i < matchList.size(); i++) {
			Integer selection = (Integer) matchList.elementAt(i);
			JyMolUtilities.label(jymol, selection.toString(), "");
			jymol.cmd.delete(selection.toString());
		}
		fireJymolClickedEvent(new JymolClickedEvent(this));
		matchList.clear();
	}

    public void setHBondVisible(boolean visible){
        showHBond = visible;
        updateHBond();
    }

    private void updateHBond() {
        jymol.cmd.delete("hbond*");
        if(receptorObj!=null&&currentHbondMolObj!=null){
            if(showHBond) {
                jymol.cmd.distance(currentHbondMolObj.getHBondName(), receptorObj.getName(), currentHbondMolObj.getName(), 2, 3.5f);
            }
            /*
            for(MolObject molObj:molObjects){
                if(molObj.isShowing()){
                    if(showHBond) {
                        jymol.cmd.distance(molObj.getHBondName(), receptorObj.getName(), molObj.getName(), 2, 3.5f);
                    }else{
                        jymol.cmd.delete(molObj.getHBondName());
                    }
                }else{
                    jymol.cmd.delete(molObj.getHBondName());
                }
            }
            */
        }
    }

    public void center(){
		if(molObject!=null){
			jymol.cmd.center(molObject.getName());
		}else{
			if(molObjects!=null&&molObjects.size()>0){
				jymol.cmd.center(molObjects.get(0).getName());
			}else {
				jymol.cmd.center("all");
			}
		}
	}

    public MolObject getReceptorObj() {
        return receptorObj;
    }

    public Vector<MolObject> getMolObjects() {
        return molObjects;
    }

    public void zoom(){
		if(molObject!=null){
			jymol.cmd.zoom(molObject.getName());
		}else{
			if(molObjects!=null&&molObjects.size()>0){
				jymol.cmd.zoom(molObjects.get(0).getName());
			}else {
				jymol.cmd.center("all");
			}
		}
	}

	public int getBuilding_mode() {
		return building_mode;
	}

	public void setBuilding_mode(int building_mode) {
        clearSelection();
		this.building_mode = building_mode;
		currentMode = BUILDING_MODE;
		if (buildMenuItem != null) {
			buildMenuItem.setSelected(true);
		}
		clearSelectionBuffer();
	}

	public String getSelectedBondParentName() {
		if (currentMode != TORSION_SELECTION_MODE) {
			return null;
		}
		if (isTorsionSelected()) {
			return bondStart.objectName;
		} else {
			return null;
		}
	}

	public void setMorphingFragment(MorphingFragment morphingFragment) {
		this.morphingFragment = morphingFragment;
	}

	public void setMorphingElement(MorphingElement morphingElement) {
		this.morphingElement = morphingElement;
	}

	public void setMorphing_linker_length(int morphing_linker_length) {
		this.morphing_linker_length = morphing_linker_length;
	}

	public void setMorphing_bond_type(String morphing_bond_type) {
		this.morphing_bond_type = morphing_bond_type;
	}


    public void setZoom(boolean zoom) {
        this.zoom = zoom;
    }

    public void addMolObject(MolObject molObject, boolean zoom){
		if(molObjectHash.containsKey(molObject.getName())){
			molObject.setHydrogen_option(hydrogen_option);
			molObjectHash.get(molObject.getName()).show(zoom);
			return;
		}
		molObject.setStyle(ligandStyle);
		molObjectHash.put(molObject.getName(),molObject);
		molObjects.add(molObject);
		molObject.setHydrogen_option(hydrogen_option);
		molObject.show(zoom);
        currentHbondMolObj = molObject;
        updateHBond();
		//clearSelection();
	}

	public void hideMolObject(String molName){
		if(molObjectHash.containsKey(molName)){
			molObjectHash.get(molName).undisplay();
		}
        currentHbondMolObj = null;
        updateHBond();
	}
    
    public void clearAllLabel(){
//        String selectionName = String.format("label_%s_%d", objectName, selection[0]);
        for(String selectionName:labelBuffer){
            String[] args = selectionName.split("_");
            String objectName = args[1];
            String selectedAtmIdx = args[2];
            int[] selection = new int[1];
            selection[0] = Integer.parseInt(selectedAtmIdx);
            jymol.cmd.selectList(selectionName,objectName,selection,0,"id",false);
            JyMolUtilities.label(jymol,selectionName,"");
            jymol.cmd.delete(selectionName);
        }
    }

	public void setReceptor(String receptorName, MolObject molObject){
		this.receptorName = receptorName;
		receptorObj = molObject;
		receptorObj.setStyle(JyMolUtilities.LINE_STYLE);
		receptorObj.setColor(receptorColor);
		receptorObj.setHydrogen_option(receptor_hydrogen_option);
		receptorIsShowing = true;
		setReceptorColor(receptorColor);
        updateHBond();
	}

	public void updateMolObject(MolObject molObject){
		updateMolObject(molObject,zoom);
	}

	public void updateMolObject(MolObject molObject,boolean zoom){
		String molName = molObject.getName();
		if(molObjectHash.containsKey(molName)){
			MolObject molObject1 = molObjectHash.get(molName);
            if(showHBond) {
                jymol.cmd.delete(molObject1.getHBondName());
            }
			molObject1.remove();
			molObjects.remove(molObject1);
			molObjectHash.remove(molName);
			addMolObject(molObject,zoom);
			molObject.setColor(molObject1.getColor());
		}else {
			addMolObject(molObject,zoom);
		}
	}

	private void removeMolObject(String molName) {
		if(molObjectHash.containsKey(molName)){
			MolObject molObject1 = molObjectHash.get(molName);
            if(showHBond) {
                jymol.cmd.delete(molObject1.getHBondName());
            }
			molObject1.remove();
			molObjects.remove(molObject1);
			molObjectHash.remove(molName);
		}
	}

	public boolean showMolObject(String molName){
		if(molObjectHash.containsKey(molName)){
            MolObject molObj = molObjectHash.get(molName);
            molObj.show(zoom);
            currentHbondMolObj = molObj;
            updateHBond();
            return true;
		}else{
			return false;
		}
	}

	public void clearAllMolObjects(){
		molObjectHash.clear();
		for(MolObject obj:molObjects){
			obj.remove();
		}
		molObjects.clear();
        if(receptorObj!=null) {
            receptorObj.remove();
        }
		clearSelection();
        currentHbondMolObj = null;
        updateHBond();
	}

	public void setMolObject(MolObject molObject) {
		if(molObject==null){
			return;
		}
		if(molObjectHash.containsKey(molObject.getName())){
			this.molObject = molObjectHash.get(molObject.getName());
		}else{
			this.molObject = molObject;
			molObjects.add(molObject);
			molObjectHash.put(molObject.getName(),molObject);
		}
		this.molObject.setHydrogen_option(hydrogen_option);
		this.molObject = molObject;
		this.molObject.updateDisplay();
        currentHbondMolObj = molObject;
	}

	public MolObject getMolObject() {
		return molObject;
	}

	public MolObject getFragObject() {
		return fragObject;
	}

	public void setFragObject(MolObject fragObject) {
		this.fragObject = fragObject;
		if (this.fragObject != null) {
			this.fragObject.setHydrogen_option(hydrogen_option);
			this.fragObject = fragObject;
			this.fragObject.updateDisplay();
		}
	}

	public int getSelectTorsionEndIndex() {
		if (torsionEnd != null) {
			return torsionEnd.member[0] - 1;
		} else {
			return -1;
		}
	}

	public int getSelectTorsionStartIndex() {
		if (torsionStart != null) {
			return torsionStart.member[0] - 1;
		} else {
			return -1;
		}
	}

	public int getSelectBondEndIndex() {
		if (bondEnd != null) {
			return bondEnd.member[0] - 1;
		} else {
			return -1;
		}
	}

	public int getSelectBondStartIndex() {
		if (bondStart != null) {
			return bondStart.member[0] - 1;
		} else {
			return -1;
		}
	}

	public boolean isTorsionSelected() {
		if (getSelectBondStartIndex() == -1 || getSelectBondEndIndex() == -1 || getSelectTorsionStartIndex() == -1 || getSelectTorsionEndIndex() == -1) {
			return false;
		} else {
			if (!bondStart.objectName.equals(bondEnd.objectName)) {
				return false;
			}
			if (!torsionStart.objectName.equals(torsionEnd.objectName)) {
				return false;
			}
			if (!torsionStart.objectName.equals(bondStart.objectName)) {
				return false;
			}
			return true;
		}
	}

	private void updateDisplay(String molName) {
        MolObject molObject = molObjectHash.get(molName);
		if (molObject!=null&&molObject.getName().equals(molName)) {
			molObject.setHydrogen_option(hydrogen_option);
			molObject.updateDisplay();
            currentHbondMolObj = molObject;
            updateHBond();
		} else if (fragObject != null && fragObject.getName().equals(molName)) {
			fragObject.setHydrogen_option(hydrogen_option);
			fragObject.updateDisplay();
		}
		if (!matchList.isEmpty()) {
			addSelection("morph", matchList);
		}
	}

	private void updateDisplay() {
		for(MolObject molObject:molObjects){
			molObject.setHydrogen_option(hydrogen_option);
			molObject.updateDisplay();
		}
		if (fragObject != null) {
			fragObject.setHydrogen_option(hydrogen_option);
			fragObject.updateDisplay();
		}
		if(receptorObj!=null){
			receptorObj.setHydrogen_option(receptor_hydrogen_option);
			receptorObj.updateDisplay();
			setReceptorStyle(receptorStyle);
		}
	}


	private void setMolString(String molName, String newMolString, boolean fireEvent) {
        MolObject molObject = molObjectHash.get(molName);
		if (molObject!=null) {
			molObject.setMolString(newMolString, fireEvent);
            updateHBond();
		}
	}

	private String getMolString(String molName) {
        MolObject molObject = molObjectHash.get(molName);
        if(molObject!=null){
            return molObject.getMolString();
        }else{
            return null;
        }
	}

	public void setReceptorColor(String receptorColor) {
		this.receptorColor = receptorColor;
		if(receptorObj!=null) {
			receptorObj.setColor(receptorColor);
			receptorObj.updateDisplay();
			setReceptorStyle(receptorStyle);
		}
//		if(receptorColor.equals(JyMolUtilities.COLOR_B_FACTOR)){
//			jymol.cmd.spectrum("b","rainbow_rev","(ss h or ss s or ss l) and e. c");
//			setReceptorStyle(receptorStyle);
//		}else if(receptorColor.equals(JyMolUtilities.COLOR_SECONDARY_STRUCTURE)){
//			jymol.cmd.color("red","ss h and e. c");
//			jymol.cmd.color("yellow","ss s and e. c");
//			jymol.cmd.color("green","ss l and e. c");
//			setReceptorStyle(receptorStyle);
//		}else{
//			if(receptorObj!=null) {
//				receptorObj.setColor(receptorColor);
//				receptorObj.updateDisplay();
//				setReceptorStyle(receptorStyle);
//			}
//		}
	}



	public void setReceptorStyle(String style){
		receptorStyle = style;
        if(style.equals(JyMolUtilities.LINE_STYLE)){
            jymol.cmd.hide(JyMolUtilities.CARTOON_STYLE, "receptor");
            jymol.cmd.hide(JyMolUtilities.SURFACE_STYLE, "receptor");
            jymol.cmd.hide(JyMolUtilities.MESH_STYLE, "receptor");
        }else if (style.equals(JyMolUtilities.PUTTY_STYLE)){
			jymol.cmd.cartoon("putty");
			jymol.cmd.hide(JyMolUtilities.SURFACE_STYLE, "receptor");
			jymol.cmd.hide(JyMolUtilities.MESH_STYLE, "receptor");
			jymol.cmd.show(JyMolUtilities.CARTOON_STYLE, "receptor");
		}else if(style.equals(JyMolUtilities.CARTOON_STYLE)){
			jymol.cmd.cartoon("automatic");
			jymol.cmd.hide(JyMolUtilities.SURFACE_STYLE, "receptor");
            jymol.cmd.hide(JyMolUtilities.MESH_STYLE, "receptor");
            jymol.cmd.show(JyMolUtilities.CARTOON_STYLE, "receptor");
        }else if(style.equals(JyMolUtilities.SURFACE_STYLE)){
			receptorObj.setColor(receptorColor);
            jymol.cmd.hide(JyMolUtilities.CARTOON_STYLE, "receptor");
            jymol.cmd.hide(JyMolUtilities.MESH_STYLE, "receptor");
            jymol.cmd.show(JyMolUtilities.SURFACE_STYLE, "receptor");
        }else if(style.equals(JyMolUtilities.MESH_STYLE)){
			receptorObj.setColor(receptorColor);
            jymol.cmd.hide(JyMolUtilities.CARTOON_STYLE, "receptor");
            jymol.cmd.hide(JyMolUtilities.SURFACE_STYLE, "receptor");
            jymol.cmd.show(JyMolUtilities.MESH_STYLE, "receptor");
        }else{
            JOptionPane.showMessageDialog(jymol,"Unsupported receptor style.");
        }
    }

	public ModelingClickConsumer(final MolViewer3D molViewer3D) {
        this.molViewer3D = molViewer3D;
		jymol = this.molViewer3D.jymol;
		molObjects = new Vector<MolObject>();
		molObjectHash = new HashMap<String, MolObject>();
		selectedMoleculeMenu = new JMenu("Selected Molecule->");
		selectedMoleculeMenu.setEnabled(false);
		JMenuItem changeColorItem = new JMenuItem("Change Color");
		molBuildingDialog = new MolBuildingDialog(this);

		changeColorItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(molObject!=null){
					String color = (String) JOptionPane.showInputDialog(molViewer3D,
							"Pick a color for selected molecule:",
							"Color Selection",
							JOptionPane.QUESTION_MESSAGE,
							null,
							JyMolUtilities.colors,
							JyMolUtilities.colors[0]);
					molObject.setColor(color);
				}
			}
		});
		JMenuItem copyAsSdfItem = new JMenuItem("Copy as SDF");
		copyAsSdfItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(molObject!=null){
					String sdfString = molObject.getMolString();
					StringSelection selection = new StringSelection(sdfString);
					Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
					clipboard.setContents(selection, selection);
					JOptionPane.showMessageDialog(molViewer3D,"Copied to clipboard.");
				}
			}
		});

        JMenuItem minimizeReceptorItem = new JMenuItem("Minimize within receptor");
        minimizeReceptorItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molObject!=null&&receptorObj!=null){
                    String sdfString = molObject.getMolString();
                    OEGraphMol oemol = OEChemFunc.getInstance().convertString2OEMol(sdfString, OEFormat.MDL);
                    //
					//todo: add minimize complex function
					//
                    OEGraphMol oemol_new = OEChemFunc.getInstance().minimizeOEMol(oemol, null/*, receptorObj.getMolString()*/);
                    String sdfString_new = OEChemFunc.getInstance().getStringFromOEMol(oemol_new);
                    molObject.setMolString(sdfString_new,true);
                    updateHBond();
                }
            }
        });


        JMenuItem minimizeItem = new JMenuItem("Minimize");
        minimizeItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molObject!=null){
                    String sdfString = molObject.getMolString();
                    OEGraphMol oemol = OEChemFunc.getInstance().convertString2OEMol(sdfString, OEFormat.MDL);
					Protonator.getInstance().protonate(oemol);
                    OEGraphMol oemol_new = OEChemFunc.getInstance().minimizeOEMol(oemol,null);
                    String sdfString_new = OEChemFunc.getInstance().getStringFromOEMol(oemol_new);
                    molObject.setMolString(sdfString_new,true);
                    updateHBond();
                }
            }
        });

        JMenuItem constrainedMinimizeItem = new JMenuItem("Minimize with Constraint ... ");
        constrainedMinimizeItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
		        if(superimposeDialog==null){
		            superimposeDialog = new CompoundInputDialog();
                }
                PropertyMolecule molecule = new PropertyMolecule(ChemFunc.getMolFromMolString(molObject.molString, OEFormat.MDL));
                superimposeDialog.setMolecule(molecule);
		        superimposeDialog.setLocationRelativeTo(molViewer3D);
		        superimposeDialog.setVisible(true);
		        if(superimposeDialog.isCommitted()){
                    Molecule mol = superimposeDialog.getMolecule();
                    if (mol != null && !mol.isEmpty()) {
                        final OEGraphMol submol = OEChemFunc.getInstance().convertChemAxonMol(mol);
                        OEGraphMol mol3d = molecule.getMol3d();
                        HashMap<Integer, Integer> matchList = OEChemFunc.getInstance().getFirstMCS(submol, mol3d, OEMCSType.Approximate);
                        if (matchList != null) {
                        	Vector<Integer> idList = new Vector<>();
                            for (OEAtomBaseIter iter = mol3d.GetAtoms(); iter.hasNext(); ) {
                                OEAtomBase atom = iter.next();
                                if (matchList.keySet().contains(atom.GetIdx())) {
                                    Integer idx = matchList.get(atom.GetIdx());
                                    idList.add(idx);
//                                    if (pred == null) {
//                                        pred = new OEHasAtomIdx(idx);
//                                    } else {
//                                        pred = new OEOrAtom(pred, new OEHasAtomIdx(idx));
//                                    }
                                }
                            }
                            OEGraphMol oeGraphMol;
                            if(receptorObj!=null) {
                                oeGraphMol = OEChemFunc.getInstance().minimizeOEMol(mol3d, idList /*, receptorObj.getMolString()*/);
                            }else{
                                oeGraphMol = OEChemFunc.getInstance().minimizeOEMol(mol3d, idList);
                            }
                            if(oeGraphMol!=null) {
                                molObject.setMolString(ChemFunc.getMolString(oeGraphMol), true);
                                updateHBond();
                            }
                        }

                    }
                }
			}
		});

		JMenuItem minimizeOPLSItem = new JMenuItem("Minimize with OPLS");
		minimizeOPLSItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(molObject!=null){
					try {
						String sdfString = molObject.getMolString();
						String minimizedMolStr = ChemFunc.minimize_mol(sdfString);
						molObject.setMolString(minimizedMolStr,true);
						updateHBond();
					} catch (MalformedURLException e1) {
						e1.printStackTrace();
					} catch (XmlRpcException e1) {
						e1.printStackTrace();
					}
				}
			}
		});

        JMenuItem minInPlaceItem = new JMenuItem("Minimize&Score (Glide)");
        minInPlaceItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(molObject!=null&&receptorObj!=null){
                    final Vector<PropertyMolecule> mols = new Vector<>();
                    PropertyMolecule propertyMolecule = new PropertyMolecule(OEChemFunc.getInstance().getOEMolFromString(molObject.getMolString()));
                    mols.add(propertyMolecule);
                    final DesignProgressMonitor progressMonitor = new DesignProgressMonitor(molViewer3D,"Progress","Progress",0,100);
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws ModelingException {
                            Vector<PropertyMolecule> v =  ChemFunc.dock(receptorName,mols,1,ChemFunc.DOCKING_METHOD_MINIMIZE,
                                    ChemFunc.DOCKING_MODE_XP,false,null,null,null,false,null);
                            return v.get(0);
                        }

                        @Override
                        protected void done() {
                            try {
                                PropertyMolecule p = (PropertyMolecule) get();
                                if(p!=null){
                                    molObject.setGlide_score(p.getProperty("Docking Score").getValue());
                                    molObject.setMolString(ChemFunc.getMolString(p.getMol3d()),true);
                                    updateHBond();
                                }
                                progressMonitor.close();
                            } catch (InterruptedException|ExecutionException e1) {
                                progressMonitor.close();
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(molViewer3D,e1.getMessage());
                            }
                        }
                    };
                    sw.execute();
                }

            }
        });


		JMenuItem rescoreItem = new JMenuItem("Rescore (Glide)");
		rescoreItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
                if(molObject!=null&&receptorObj!=null){
                    final Vector<PropertyMolecule> mols = new Vector<>();
                    PropertyMolecule propertyMolecule = new PropertyMolecule(OEChemFunc.getInstance().getOEMolFromString(molObject.getMolString()));
                    mols.add(propertyMolecule);
                    final DesignProgressMonitor progressMonitor = new DesignProgressMonitor(molViewer3D,"Progress","Progress",0,100);
                    progressMonitor.setProgress(DesignProgressMonitor.INDETERMINATE);
                    SwingWorker sw = new SwingWorker() {
                        @Override
                        protected Object doInBackground() throws Exception {
                            Vector<PropertyMolecule> v =  ChemFunc.dock(receptorName,mols,1,ChemFunc.DOCKING_METHOD_SCOREONLY,
                                    ChemFunc.DOCKING_MODE_XP,false,null,null,null,false,null);
                            return v.get(0);
                        }

                        @Override
                        protected void done() {
                            try {
                                PropertyMolecule p = (PropertyMolecule) get();
                                if(p!=null&&molObject!=null){
                                    molObject.setGlide_score(p.getProperty("Docking Score").getValue());
                                    molObject.setMolString(ChemFunc.getMolString(p.getMol3d()),true);
                                    updateHBond();
                                }
                                progressMonitor.close();
                            } catch (InterruptedException | ExecutionException e1) {
                                progressMonitor.close();
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(molViewer3D,e1.getMessage());
                            }
                        }
                    };
                    sw.execute();
                }
			}
		});


		selectedMoleculeMenu.add(changeColorItem);
		selectedMoleculeMenu.add(copyAsSdfItem);
        selectedMoleculeMenu.add(minimizeItem);
        selectedMoleculeMenu.add(constrainedMinimizeItem);
        selectedMoleculeMenu.add(minimizeOPLSItem);
        selectedMoleculeMenu.add(minimizeReceptorItem);
        selectedMoleculeMenu.add(rescoreItem);
        selectedMoleculeMenu.add(minInPlaceItem);

		ButtonGroup group = new ButtonGroup();
		buildMenuItem = new JRadioButtonMenuItem("Building Mode");
		buildMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
				if (menuItem.isSelected()) {
					currentMode = BUILDING_MODE;
					clearSelectionBuffer();
					clearSelection();
					molBuildingDialog.setVisible(true);
				}
			}
		});


        JRadioButtonMenuItem labelMenuItem = new JRadioButtonMenuItem("Label Mode");
        labelMenuItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
                if (menuItem.isSelected()) {
                    currentMode = LABELING_MODE;
                    clearSelectionBuffer();
                    clearSelection();
                }
            }
        });

		JRadioButtonMenuItem distanceMenuItem = new JRadioButtonMenuItem("Distance Mode");
		distanceMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
				if (menuItem.isSelected()) {
					currentMode = DISTANCE_MODE;
					clearSelectionBuffer();
					clearSelection();
				}
			}
		});
		JRadioButtonMenuItem angleMenuItem = new JRadioButtonMenuItem("Angle Mode");
		angleMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
				if (menuItem.isSelected()) {
					currentMode = ANGLE_MODE;
					clearSelectionBuffer();
					clearSelection();
				}
			}
		});
		JRadioButtonMenuItem dihedralsMenuItem = new JRadioButtonMenuItem("Dihedrals Mode");
		dihedralsMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
				if (menuItem.isSelected()) {
					currentMode = DIHEDRAL_MODE;
					clearSelectionBuffer();
					clearSelection();
				}
			}
		});
		final JRadioButtonMenuItem selectionMenuItem = new JRadioButtonMenuItem("Selection Mode");
		selectionMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
				if (menuItem.isSelected()) {
					currentMode = SELECTION_MODE;
					clearSelectionBuffer();
					clearSelection();
				}
			}
		});
		selectionMenuItem.setSelected(true);
		molBuildingDialog.addWindowListener(new WindowListener() {
			@Override
			public void windowOpened(WindowEvent e) {

			}

			@Override
			public void windowClosing(WindowEvent e) {
				clearSelectionBuffer();
				clearSelection();
				buildMenuItem.setSelected(false);
				selectionMenuItem.setSelected(true);
				currentMode = SELECTION_MODE;
			}

			@Override
			public void windowClosed(WindowEvent e) {
			}

			@Override
			public void windowIconified(WindowEvent e) {

			}

			@Override
			public void windowDeiconified(WindowEvent e) {

			}

			@Override
			public void windowActivated(WindowEvent e) {

			}

			@Override
			public void windowDeactivated(WindowEvent e) {

			}
		});


		JRadioButtonMenuItem torsionSelectionMenuItem = new JRadioButtonMenuItem("Torsion Selection Mode");
		torsionSelectionMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JRadioButtonMenuItem menuItem = (JRadioButtonMenuItem) e.getSource();
				if (menuItem.isSelected()) {
					currentMode = TORSION_SELECTION_MODE;
					clearSelectionBuffer();
					clearSelection();
				}
			}
		});

		buildMenuItem.setSelected(false);

		JMenuItem clearAllMenuItem = new JMenuItem("Clear all monitors");
		clearAllMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (JOptionPane.showConfirmDialog(jymol, "Are you sure?") == JOptionPane.OK_OPTION) {
					clearAllMonitors();
				}
			}
		});

        JMenuItem clearAllLabelItem = new JMenuItem("Clear all labels");
        clearAllLabelItem.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (JOptionPane.showConfirmDialog(jymol, "Are you sure?") == JOptionPane.OK_OPTION) {
                    clearAllLabel();
                }
            }
        });

		group.add(buildMenuItem);
		group.add(selectionMenuItem);
		group.add(distanceMenuItem);
		group.add(angleMenuItem);
		group.add(dihedralsMenuItem);
		group.add(torsionSelectionMenuItem);
        group.add(labelMenuItem);


		ButtonGroup group2 = new ButtonGroup();
		final JMenu showHydrogenMenu = new JMenu("Ligand hydrogen style:");
		JRadioButtonMenuItem showNoHydro = new JRadioButtonMenuItem("Show No Hydrogen", true);
		showNoHydro.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				hydrogen_option = SHOW_NO_HYDROGEN;
				updateDisplay();
			}
		});
		JRadioButtonMenuItem showPolarHydro = new JRadioButtonMenuItem("Show Polar Hydrogen Only", false);
		showPolarHydro.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				hydrogen_option = SHOW_POLAR_HYDROGEN;
				updateDisplay();
			}
		});
		JRadioButtonMenuItem showAllHydro = new JRadioButtonMenuItem("Show All Hydrogen", false);
		showAllHydro.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				hydrogen_option = SHOW_ALL_HYDROGEN;
				updateDisplay();
			}
		});
		group2.add(showNoHydro);
		group2.add(showPolarHydro);
		group2.add(showAllHydro);
		showHydrogenMenu.add(showNoHydro);
		showHydrogenMenu.add(showPolarHydro);
		showHydrogenMenu.add(showAllHydro);

		ButtonGroup group4 = new ButtonGroup();
		final JMenu showReceptorHydrogenMenu = new JMenu("Receptor hydrogen style:");
		JRadioButtonMenuItem showNoReceptorHydrogen = new JRadioButtonMenuItem("Show No Hydrogen",true);
		showNoReceptorHydrogen.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				receptor_hydrogen_option = SHOW_NO_HYDROGEN;
				updateDisplay();
			}
		});
		JRadioButtonMenuItem showPolarReceptorHydrogen = new JRadioButtonMenuItem("Show Polar Hydrogen",false);
		showPolarReceptorHydrogen.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				receptor_hydrogen_option = SHOW_POLAR_HYDROGEN;
				updateDisplay();
			}
		});
		JRadioButtonMenuItem showAllReceptorHydrogen = new JRadioButtonMenuItem("Show All Hydrogen",false);
		showAllReceptorHydrogen.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				receptor_hydrogen_option = SHOW_ALL_HYDROGEN;
				updateDisplay();
			}
		});
		group4.add(showNoReceptorHydrogen);
		group4.add(showPolarReceptorHydrogen);
		group4.add(showAllReceptorHydrogen);
		showReceptorHydrogenMenu.add(showNoReceptorHydrogen);
		showReceptorHydrogenMenu.add(showPolarReceptorHydrogen);
		showReceptorHydrogenMenu.add(showAllReceptorHydrogen);


		menu.add(buildMenuItem);
		menu.add(selectedMoleculeMenu);
		menu.add(showHydrogenMenu);
		menu.add(showReceptorHydrogenMenu);
		menu.add(selectionMenuItem);
		menu.add(distanceMenuItem);

		JMenu monitorMenu = new JMenu("Other Picking Modes ...");
		monitorMenu.add(angleMenuItem);
		monitorMenu.add(dihedralsMenuItem);
		monitorMenu.add(torsionSelectionMenuItem);
        monitorMenu.add(labelMenuItem);

		menu.add(monitorMenu);
		menu.add(clearAllMenuItem);
        menu.add(clearAllLabelItem);

		JMenuItem bgcolor = new JMenuItem("Change Backgroud Color ...");
		bgcolor.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				Color bgColor = JColorChooser.showDialog(jymol, "Choose Background Color", Color.black);
				if (bgColor != null) {
					float[] rgbf = bgColor.getRGBColorComponents(null);
					jymol.cmd.set("bg_rgb", JyMolUtilities.convertColorRGB(rgbf));
				}
			}
		});
		menu.add(bgcolor);

		ButtonGroup buttonGroup3 = new ButtonGroup();
		final JMenu receptorStyleMenu = new JMenu("Receptor Styles:");
		JRadioButtonMenuItem lineStyleItem = new JRadioButtonMenuItem("Line");
		lineStyleItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
                setReceptorStyle(JyMolUtilities.LINE_STYLE);
			}
		});

		JRadioButtonMenuItem meshStyleItem = new JRadioButtonMenuItem("Mesh");
		meshStyleItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
                setReceptorStyle(JyMolUtilities.MESH_STYLE);
			}
		});
		JRadioButtonMenuItem cartoonStyleItem = new JRadioButtonMenuItem("Cartoon");
		cartoonStyleItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
                setReceptorStyle(JyMolUtilities.CARTOON_STYLE);
			}
		});
		JRadioButtonMenuItem surfaceStyleItem = new JRadioButtonMenuItem("Surface");
		surfaceStyleItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
                setReceptorStyle(JyMolUtilities.SURFACE_STYLE);
			}
		});
		buttonGroup3.add(lineStyleItem);
		buttonGroup3.add(meshStyleItem);
		buttonGroup3.add(cartoonStyleItem);
		buttonGroup3.add(surfaceStyleItem);
		receptorStyleMenu.add(lineStyleItem);
		receptorStyleMenu.add(meshStyleItem);
		receptorStyleMenu.add(cartoonStyleItem);
		receptorStyleMenu.add(surfaceStyleItem);
		menu.add(receptorStyleMenu);


		JMenuItem rayTracingMenuItem = new JMenuItem("Generate screenshot as is ...");
		rayTracingMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent actionEvent) {
				if (jymol != null) {
					jymol.cmd.ray();
				}
			}
		});
		menu.add(new JSeparator());
		menu.add(rayTracingMenuItem);
		JMenuItem powerpointMenuItem = new JMenuItem("Create a picture for powerpoint ...");
		powerpointMenuItem.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				float[] rgbf = Color.WHITE.getRGBColorComponents(null);
				jymol.cmd.set("bg_rgb", JyMolUtilities.convertColorRGB(rgbf));
                String[] objectList = jymol.cmd.getObjectNames();
                if(objectList!=null&&objectList.length>0){
                    List<String> names = Arrays.asList(objectList);
                    if(names.contains("receptor")){
                        jymol.cmd.show("sticks", "all, not receptor");
                    }else{
                        jymol.cmd.show("sticks", "all");
                    }
                }
				/*
								if (receptorIsShowing) {
									jymol.cmd.show("sticks", "all, not receptor");
								} else {
									jymol.cmd.show("sticks", "all");
								}
								*/
				jymol.cmd.ray();
			}
		});
		menu.add(powerpointMenuItem);

		final JCheckBoxMenuItem stereoMenuItem = new JCheckBoxMenuItem("Stereo Mode (Red-blue)", ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("3D-Glasses-icon.png"))));
		stereoMenuItem.setState(false);
		stereoMenuItem.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if(stereoMenuItem.isSelected()){
					jymol.cmd.set("stereo_mode","10");
					jymol.cmd.set("stereo","1");
				}else{
					jymol.cmd.set("stereo","0");
				}
			}
		});
		menu.add(stereoMenuItem);

		updateDisplay();
		addJymolClickedEventListener(this);
	}

	public JPopupMenu getPopupMenu() {
		return menu;
	}

	public void clearSelectionBuffer() {
		for (int i = 0; i < selectionBuffer.size(); i++) {
			JymolSelection jymolSelection = selectionBuffer.elementAt(i);
			jymolSelection.remove();
		}
		selectionBuffer.clear();
	}

	private void clearSelectedTorsion() {
		if (bondEnd != null) {
			bondEnd.remove();
			bondEnd = null;
		}
		if (bondStart != null) {
			bondStart.remove();
			bondStart = null;
		}
		if (torsionEnd != null) {
			torsionEnd.remove();
			torsionEnd = null;
		}
		if (torsionStart != null) {
			torsionStart.remove();
			torsionStart = null;
		}
	}

	private String getNewName(String s) {
		return s;
	}


	public void processClick(HashMap hashMap) {
		if (hashMap == null || jymol == null) return;
		String clickInfo = (String) hashMap.get("click");
		String modKeys = (String) hashMap.get("mod_keys");
		if (clickInfo.equals("single_right") || (clickInfo.equals("single_left") && modKeys.equals("ctrl"))) {
			int x = Integer.parseInt((String) hashMap.get("x"));
			int y = Integer.parseInt((String) hashMap.get("y"));
			menu.show(jymol, x, y);
		} else if (clickInfo.equals("single_left")) {
			String objectType = (String) hashMap.get("type");
			if (objectType == null){
                clearSelection();
                return;
            }
			Object idObj = hashMap.get("id");
			if(currentMode == LABELING_MODE){
				if(idObj==null){
					return;
				}
//                type object:molecule
//                click single_left
//                mod_keys
//                name N
//                x 795
//                rank 4708
//                y 396
//                id 4710
//                resn THR
//                state 1
//                resi 94
//                object receptor
                String atomName = (String) hashMap.get("name");
                String residue_type = (String)hashMap.get("resn");
                String chain_name = (String) hashMap.get("chain");
                String residue_number = (String) hashMap.get("resi");
                String objectName = (String)hashMap.get("object");
                String segment = (String)hashMap.get("segi");
                int[] selection = new int[1];
                selection[0] = Integer.parseInt((String) idObj);
                String selectionName = String.format("label_%s_%d", objectName, selection[0]);
                if(!labelBuffer.contains(selectionName)){
                    labelBuffer.add(selectionName);
                    jymol.cmd.selectList(selectionName,objectName,selection,0,"id",false);
                    JyMolUtilities.label(jymol,selectionName,String.format("%s%s %s",residue_type,residue_number,atomName));
                }else{
                    labelBuffer.remove(selectionName);
                    jymol.cmd.selectList(selectionName,objectName,selection,0,"id",false);
                    JyMolUtilities.label(jymol,selectionName,"");
                }
                jymol.cmd.delete(selectionName);
                return;
            }

			if (currentMode == SELECTION_MODE) {
				if(objectType.equals("object:molecule")){
					String molName = (String) hashMap.get("object");
					if(molObjectHash.containsKey(molName)) {
						String molString = molObjectHash.get(molName).getMolString();
						OEGraphMol mol = OEChemFunc.getInstance().getOEMolFromString(molString);
						if (mol != null) {
							Vector<Integer> selection = new Vector<Integer>();
							for (OEAtomBase atm : mol.GetAtoms()) {
								selection.add(atm.GetIdx());
							}
							addSelection(molName, selection);
						}
					}
				}else{
                    clearSelection();
                }
				return;
			}

			if (!objectType.equals("object:molecule")) {
				return;
			}

			if(idObj==null){
				return;
			}

			String molName = (String) hashMap.get("object");
			Integer idx = Integer.valueOf((String) idObj);
			JymolSelection selection = new JymolSelection(molName, idx);
			if (selectionBuffer.contains(selection)) {
				selectionBuffer.remove(selection);
				selection.remove();
			} else {
				selectionBuffer.add(selection);
			}
			/*
						System.out.println("___________________________________________________");
						for (int i = 0; i < selectionBuffer.size(); i++) {
							JymolSelection jymolSelection = selectionBuffer.elementAt(i);
							System.out.println("JYMOL SELECTION:" + jymolSelection.toString());
						}
						System.out.println("___________________________________________________");
						*/
			switch (currentMode) {
				case BUILDING_MODE:
					String newMolString = null;
					if(!molObjectHash.containsKey(molName)){
						return;
					}
                    try {
						MolObject oldObj = molObjectHash.get(molName);
						undoStack.push(new MolObject(oldObj));
						switch (building_mode) {
							case BUILDING_SUBMODE_ATTACH:
								newMolString = MoleculeMorpher.attachGroupToMolecule(getMolString(molName), idx, morphingFragment);
								setMolString(molName, newMolString, true);
								clearSelectionBuffer();
								break;
							case BUILDING_SUBMODE_DELETE:
								newMolString = MoleculeMorpher.deleteAtomFromMolecule(getMolString(molName), idx);
								setMolString(molName, newMolString, true);
								break;
							case BUILDING_SUBMODE_LINK:
								if (selectionBuffer.size() == 2) {
									int idx1 = selectionBuffer.get(0).getMember();
									int idx2 = selectionBuffer.get(1).getMember();
									if (selectionBuffer.get(0).getObjectName().equals(selectionBuffer.get(1).getObjectName())) {
										newMolString = MoleculeMorpher.linkAtoms(getMolString(molName), idx1, idx2, morphing_linker_length);
										setMolString(molName, newMolString, true);
									} else {
										String parentString = getMolString(molObject.getName());
										String fragmentString = getMolString(fragObject.getName());
										String name1 = selectionBuffer.get(0).getObjectName();
										String name2 = selectionBuffer.get(1).getObjectName();
										if (name1.equals(molObject.getName()) && name2.equals(fragObject.getName())) {
											newMolString = MoleculeMorpher.linkMolecule(parentString, fragmentString, idx1, idx2, morphing_linker_length);
										} else if (name1.equals(fragObject.getName()) && name2.equals(molObject.getName())) {
											newMolString = MoleculeMorpher.linkMolecule(parentString, fragmentString, idx2, idx1, morphing_linker_length);
										}
										if (newMolString != null) {
											setMolString(molObject.getName(), newMolString, true);
										}else{
											undoStack.pop();
										}
									}
									clearSelectionBuffer();
								}
								break;
							case BUILDING_SUBMODE_MODIFY:
								newMolString = MoleculeMorpher.modifyAtom(getMolString(molName), idx, morphingElement);
								setMolString(molName, newMolString, true);
								clearSelectionBuffer();
								break;
							case BUILDING_SUBMODE_MODIFY_BOND:
								if (selectionBuffer.size() == 2 && selectionBuffer.get(0).getObjectName().equals(selectionBuffer.get(1).getObjectName())) {
									int idx1 = selectionBuffer.get(0).getMember();
									int idx2 = selectionBuffer.get(1).getMember();
									newMolString = MoleculeMorpher.modifyBond(getMolString(molName), idx1, idx2, morphing_bond_type);
									setMolString(molName, newMolString, true);
									clearSelectionBuffer();
								}else{
									undoStack.pop();
								}
								break;
							case BUILDING_SUBMODE_DELETE_BOND:
								if (selectionBuffer.size() == 2) {
									if (selectionBuffer.get(0).getObjectName().equals(selectionBuffer.get(1).getObjectName())) {
										int idx1 = selectionBuffer.get(0).getMember();
										int idx2 = selectionBuffer.get(1).getMember();
										newMolString = MoleculeMorpher.modifyBond(getMolString(molName), idx1, idx2, "none");
										setMolString(molName, newMolString, true);
									}else{
										undoStack.pop();
									}
									clearSelectionBuffer();
								}
								break;
							case BUILDING_SUBMODE_ROTATE_BOND:
								if (selectionBuffer.size() == 4) {
									addTorsionSelection(selectionBuffer.get(0), selectionBuffer.get(1), selectionBuffer.get(2), selectionBuffer.get(3));
									//clearSelectionBuffer();
								}
								break;
							default:
								undoStack.pop();
								break;
						}
					} catch (final Exception e) {
						undoStack.pop();
						SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(jymol.getParent(), e.getMessage());
								clearSelectionBuffer();
							}
						});
					}
					break;
				case DISTANCE_MODE:
					if (selectionBuffer.size() == 2) {
						String[] pointSelections = new String[2];
						pointSelections[0] = selectionBuffer.get(0).selectionName;
						pointSelections[1] = selectionBuffer.get(1).selectionName;
						addDistance(pointSelections);
						clearSelectionBuffer();
					}
					break;
				case ANGLE_MODE:
					if (selectionBuffer.size() == 3) {
						String[] pointSelections = new String[3];
						pointSelections[0] = selectionBuffer.get(0).selectionName;
						pointSelections[1] = selectionBuffer.get(1).selectionName;
						pointSelections[2] = selectionBuffer.get(2).selectionName;
						addAngle(pointSelections);
						clearSelectionBuffer();
					}
					break;
				case DIHEDRAL_MODE:
					if (selectionBuffer.size() == 4) {
						String[] pointSelections = new String[4];
						pointSelections[0] = selectionBuffer.get(0).selectionName;
						pointSelections[1] = selectionBuffer.get(1).selectionName;
						pointSelections[2] = selectionBuffer.get(2).selectionName;
						pointSelections[3] = selectionBuffer.get(3).selectionName;
						addDihedral(pointSelections);
						clearSelectionBuffer();
					}
					break;
				case TORSION_SELECTION_MODE:

					/*
										if(selectionBuffer!=null&&selectionBuffer.size()==1){
											if(selectionBuffer.get(0).equals(bondStart)||selectionBuffer.get(0).equals(bondEnd)){
												if(bondStart!=null){
													bondStart.remove();
													bondStart = null;
												}
												if(bondEnd!=null){
													bondEnd.remove();
													bondEnd = null;
												}
												clearSelectionBuffer();
												return;
											}
										}
										*/
					if (selectionBuffer.size() == 4) {
						addTorsionSelection(selectionBuffer.get(0), selectionBuffer.get(1), selectionBuffer.get(2), selectionBuffer.get(3));
						selectionBuffer.clear();
					}
					break;
				default:
					break;
			}
		}
		return;
	}

	public void clearAllMonitors() {
		numDistances = 0;
		numAngles = 0;
		numDihedrals = 0;

		for (Iterator<String> iterator = psName2dName.values().iterator(); iterator.hasNext(); ) {
			String monitorName = iterator.next();
			jymol.cmd.delete(monitorName);
		}
		psName2dName.clear();
	}

	void addTorsionSelection(JymolSelection torsionStart1, JymolSelection bondStart1, JymolSelection bondEnd1, JymolSelection torsionEnd1) {
		if (torsionStart1 == null || torsionEnd1 == null || bondStart1 == null || bondEnd1 == null) {
			return;
		}
		if (!torsionStart1.getObjectName().equals(torsionEnd1.getObjectName()) || (!torsionEnd1.getObjectName().equals(bondStart1.getObjectName())) || (!bondStart1.getObjectName().equals(bondEnd1.getObjectName()))) {
			return;
		}
		if (torsionStart != null && !torsionStart.equals(torsionStart1)) {
			torsionStart.remove();
		}
		if (bondStart != null && !bondStart.equals(bondStart1)) {
			bondStart.remove();
		}
		if (bondEnd != null && !bondEnd.equals(bondEnd1)) {
			bondEnd.remove();
		}
		if (torsionEnd != null && !torsionEnd.equals(torsionEnd1)) {
			torsionEnd.remove();
		}
		torsionStart = new JymolSelection(torsionStart1.getObjectName(), torsionStart1.getMember());
		JyMolUtilities.label(jymol, torsionStart.selectionName, "A");
		bondStart = new JymolSelection(bondStart1.getObjectName(), bondStart1.getMember());
		JyMolUtilities.label(jymol, bondStart.selectionName, "B");
		bondEnd = new JymolSelection(bondEnd1.getObjectName(), bondEnd1.getMember());
		JyMolUtilities.label(jymol, bondEnd.selectionName, "C");
		torsionEnd = new JymolSelection(torsionEnd1.getObjectName(), torsionEnd1.getMember());
		JyMolUtilities.label(jymol, torsionEnd.selectionName, "D");

	}

	void addDistance(String[] pointSelections) {
		if (!pointSelections[0].equals(pointSelections[1])) {
			boolean deleted = false;
			for (int i = 0; i < 2 && !deleted; i++) {
				for (int j = 0; j < 2 && !deleted; j++) {
					if (i != j) {
						String name = getNewName(pointSelections[i] + "___" + pointSelections[j]);
						if (psName2dName.containsKey(name)) {
							jymol.cmd.delete(psName2dName.get(name));
							psName2dName.remove(name);
							deleted = true;
							numDistances--;
						}
					}
				}
			}
			if (!deleted) {
				jymol.cmd.distance("d_" + numDistances, pointSelections[0], pointSelections[1], 0, -1.0f, true, false, 0, 0, false);
				psName2dName.put(getNewName(pointSelections[0] + "___" + pointSelections[1]), "d_" + numDistances);
				numDistances++;
			}

		}
	}

	void addAngle(String[] pointSelections) {
		if ((!pointSelections[0].equals(pointSelections[1])) &&
				(!pointSelections[1].equals(pointSelections[2])) &&
				(!pointSelections[0].equals(pointSelections[2]))) {
			String name1 = getNewName(pointSelections[0] + "___" + pointSelections[1] + "___" + pointSelections[2]);
			String name2 = getNewName(pointSelections[2] + "___" + pointSelections[1] + "___" + pointSelections[0]);
			if (psName2dName.containsKey(name1)) {
				jymol.cmd.delete(psName2dName.get(name1));
				psName2dName.remove(name1);
				numAngles--;
			} else if (psName2dName.containsKey(name2)) {
				jymol.cmd.delete(psName2dName.get(name2));
				psName2dName.remove(name2);
				numAngles--;
			} else {
				jymol.cmd.angle("a_" + numAngles, pointSelections[0], pointSelections[1], pointSelections[2], 0, true, false, 0, 0, false);
				System.out.println("iname is " + "a_" + numAngles);
				psName2dName.put(getNewName(pointSelections[0] + "___" + pointSelections[1] + "___" + pointSelections[2]), "a_" + numAngles);
				numAngles++;
			}
		}
	}

	void addDihedral(String[] pointSelections) {
		if ((!pointSelections[0].equals(pointSelections[1])) &&
				(!pointSelections[1].equals(pointSelections[2])) &&
				(!pointSelections[0].equals(pointSelections[2])) &&
				(!pointSelections[0].equals(pointSelections[3])) &&
				(!pointSelections[1].equals(pointSelections[3])) &&
				(!pointSelections[2].equals(pointSelections[3]))) {
			String name1 = getNewName(pointSelections[0] + "___" + pointSelections[1] + "___" + pointSelections[2] + "___" + pointSelections[3]);
			String name2 = getNewName(pointSelections[3] + "___" + pointSelections[2] + "___" + pointSelections[1] + "___" + pointSelections[0]);
			if (psName2dName.containsKey(name1)) {
				jymol.cmd.delete(psName2dName.get(name1));
				psName2dName.remove(name1);
				numDihedrals--;
			} else if (psName2dName.containsKey(name2)) {
				jymol.cmd.delete(psName2dName.get(name2));
				psName2dName.remove(name2);
				numDihedrals--;
			} else {
				jymol.cmd.dihedral("dh_" + numDihedrals, pointSelections[0], pointSelections[1], pointSelections[2], pointSelections[3], 0, true, false, 0, 0, false);
				System.out.println("iname is " + "dh_" + numDihedrals);
				System.out.println(getNewName(pointSelections[0] + "___" + pointSelections[1] + "___" + pointSelections[2] + "___" + pointSelections[3]));
				psName2dName.put(getNewName(pointSelections[0] + "___" + pointSelections[1] + "___" + pointSelections[2] + "___" + pointSelections[3]), "dh_" + numDihedrals);
				numDihedrals++;
			}
		}
	}

	public void rotateMolecule(double value, boolean fireEvent) throws Exception {
		if (building_mode == BUILDING_SUBMODE_ROTATE_BOND) {

			int bgnTorsionIdx = getSelectTorsionStartIndex();
			int bgnIdx = getSelectBondStartIndex();
			int endIdx = getSelectBondEndIndex();
			int endTorsionIdx = getSelectTorsionEndIndex();
			if (bgnTorsionIdx == -1 || bgnIdx == -1 || endIdx == -1 || endTorsionIdx == -1) {
				throw new Exception("No Torsion defined.");
			}
			String newMolString = OEChemFunc.getInstance().rotateBond(getMolString(torsionStart.getObjectName()), bgnTorsionIdx, bgnIdx, endIdx, endTorsionIdx, value);
			setMolString(torsionStart.getObjectName(), newMolString, false);
			updateDisplay(torsionStart.getObjectName());
			if (fireEvent) {
				jymol.firePropertyChange("rotate", true, false);
			}
		}
	}

	public void setAllStyles(String style) {
		ligandStyle = style;
		for(MolObject obj:molObjects){
			obj.setStyle(style);
		}
	}

	@Override
	public void jymolClicked(JymolClickedEvent event) {
		if(molObject!=null){
			selectedMoleculeMenu.setEnabled(true);
		}else{
			selectedMoleculeMenu.setEnabled(false);
		}
	}


	class JymolSelection extends Object {
		String selectionName;
		String objectName;
		int[] member = new int[1];

		public JymolSelection(String name, int member_idx) {
			selectionName = name + member_idx;
			objectName = name;
			member[0] = member_idx;
			jymol.cmd.selectList(selectionName, objectName, member);
			JyMolUtilities.label(jymol, selectionName, "*");
		}

		public int getMember() {
			return member[0];
		}

		public String getObjectName() {
			return objectName;
		}

		public String getSelectionName() {
			return selectionName;
		}

		public void remove() {
			JyMolUtilities.label(jymol, selectionName, "");
			jymol.cmd.delete(selectionName);
		}

		public boolean equals(Object obj) {
			if (!(obj instanceof JymolSelection)) return false;
			if (this.selectionName.equals(((JymolSelection) obj).selectionName)) return true;
			return false;
		}

		public String toString() {
			return String.format("%s_%s_%d", selectionName, objectName, member[0]);
		}
	}


}
