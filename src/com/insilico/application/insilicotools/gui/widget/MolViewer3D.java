package com.insilico.application.insilicotools.gui.widget;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.util.JyMolUtilities;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemWebLicenseInstaller;
import com.schrodinger.jymol.JyMol;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;

import javax.swing.*;
import javax.swing.border.BevelBorder;
import javax.swing.border.LineBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.*;

/**
 * Created by jfeng1 on 5/2/16.
 */
public class MolViewer3D extends JPanel{
    JyMol jymol;
    ModelingClickConsumer clickConsumer;
    boolean useReceptor = false;
    JCheckBox zoomBx;
    JCheckBox hbondBx;
    ModelingImageConsumer imageConsumer;

    public MolViewer3D() {
        this(false,false);
    }

    public void setZoom(boolean zoom) {
        if(zoomBx!=null){
            zoomBx.setSelected(zoom);
            clickConsumer.setZoom(zoom);
        }
    }

    public MolViewer3D(boolean useBtnPanel, boolean useReceptor){
        super(new BorderLayout());
        this.useReceptor = useReceptor;
        initialize();
        if(useBtnPanel){
            add(buildBtnPanel(useReceptor),BorderLayout.SOUTH);
        }
        setBorder(new BevelBorder(BevelBorder.RAISED));
    }

    private void initialize() {
        jymol = new JyMol(JyMol.ENABLE_STEREO);
        jymol.cmd.feedback(2,0,0xff);
        imageConsumer = new ModelingImageConsumer(this);
        jymol.addImageConsumer(imageConsumer);
        JyMolUtilities.setHighQuality(jymol);
        JPanel p = new JPanel();
        p.setLayout(new GridLayout(1, 1, 3, 3));
        p.add(jymol);
        p.setMinimumSize(new Dimension(0, 0));
        p.setMaximumSize(new Dimension(0, 0));
        clickConsumer = new ModelingClickConsumer(this);
        jymol.addClickConsumer(clickConsumer);
        add(p,BorderLayout.CENTER);
    }


    public JPanel buildBtnPanel(boolean useReceptor) {
        JPanel p = new JPanel();
        final JComboBox ligStyleCB = new JComboBox(new String[]{JyMolUtilities.STICK_STYLE,JyMolUtilities.LINE_STYLE,JyMolUtilities.SURFACE_STYLE,JyMolUtilities.CPK_STYLE});
        ligStyleCB.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                clickConsumer.setAllStyles((String)ligStyleCB.getSelectedItem());
            }
        });
        JPanel p1 = new JPanel();
        p1.add(new JLabel("Ligand Style:"));
        p1.add(ligStyleCB);
        p1.setBorder(new LineBorder(Color.black,1));
        p.add(p1);

        if(useReceptor) {
            JPanel p2 = new JPanel();
            p2.add(new JLabel("Receptor Style:"));
            final JComboBox receptorStyleCB = new JComboBox(new String[]{JyMolUtilities.LINE_STYLE, JyMolUtilities.MESH_STYLE, JyMolUtilities.SURFACE_STYLE, JyMolUtilities.CARTOON_STYLE, JyMolUtilities.PUTTY_STYLE});
            receptorStyleCB.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    clickConsumer.setReceptorStyle((String)receptorStyleCB.getSelectedItem());
                }
            });
            p2.add(receptorStyleCB);

            p2.add(new JToolBar.Separator());
            Vector<String> receptorColors = new Vector<String>(Arrays.asList(JyMolUtilities.colors));
            receptorColors.add(JyMolUtilities.COLOR_B_FACTOR);
            receptorColors.add(JyMolUtilities.COLOR_SECONDARY_STRUCTURE);
            final JComboBox receptorColorCB = new JComboBox(receptorColors);
            receptorColorCB.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    clickConsumer.setReceptorColor((String)receptorColorCB.getSelectedItem());
                }
            });
            p2.add(new JLabel("Color:"));
            p2.add(receptorColorCB);
            p2.add(new JToolBar.Separator());

            hbondBx = new JCheckBox("HBond",true);
            hbondBx.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    clickConsumer.setHBondVisible(hbondBx.isSelected());
                }
            });
            p2.add(hbondBx);
            p2.setBorder(new LineBorder(Color.black, 1));
            p.add(p2);
        }

        JPanel p3 = new JPanel();
        JButton centerBtn = new JButton("Recenter");
        centerBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                clickConsumer.center();
                clickConsumer.zoom();
            }
        });
        p3.add(centerBtn);

        zoomBx = new JCheckBox("Auto Zoom");
        zoomBx.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                clickConsumer.setZoom(zoomBx.isSelected());
            }
        });
        p3.add(zoomBx);

//        JButton rezoomBtn = new JButton("Zoom");
//        rezoomBtn.addActionListener(new ActionListener() {
//            @Override
//            public void actionPerformed(ActionEvent e) {
//                clickConsumer.zoom();
//            }
//        });
//        p3.add(rezoomBtn);
        p3.setBorder(new LineBorder(Color.black,1));
        p.add(p3);
        return p;
    }

    public void clear(){
        clickConsumer.clearAllMolObjects();
        jymol.cmd.reinitialize();
        JyMolUtilities.setHighQuality(jymol);
    }

    public boolean updateLigand(PropertyMolecule mol){
        OEGraphMol m3d = mol.getMol3d();
        if (m3d != null) {
            String molString = ChemFunc.getMolString(m3d);
            String molName = mol.getUniqName();
            int color_id = new Random().nextInt(JyMolUtilities.colors.length);
            clickConsumer.updateMolObject(new MolObject(molName, molString, JyMolUtilities.colors[color_id], jymol));
            return true;
        } else {
            return false;
        }
    }

    public boolean addLigand(PropertyMolecule mol){
        if(!clickConsumer.showMolObject(mol.getUniqName())) {
            OEGraphMol m3d = mol.getMol3d();
            if (m3d != null) {
                String molString = ChemFunc.getMolString(m3d);
                String molName = mol.getUniqName();
                int color_id = new Random().nextInt(JyMolUtilities.colors.length);
                MolObject molObject = new MolObject(molName, molString, JyMolUtilities.colors[color_id], jymol);
                clickConsumer.addMolObject(molObject,clickConsumer.zoom);
                return true;
            } else {
                return false;
            }
        }
        return true;
    }

    public boolean setReceptor(String receptorName, String pdbStr){
        if(pdbStr!=null&&pdbStr.length()>0) {
            clickConsumer.setReceptor(receptorName, new MolObject("receptor", pdbStr, JyMolUtilities.colors[0], jymol, "pdb"));
            return true;
        }else{
            return false;
        }
    }



    public void hideLigand(PropertyMolecule mol){
        if(mol!=null){
            clickConsumer.hideMolObject(mol.getUniqName());
        }
    }

    public static void main(String[] args) {
        try {
            OEChemWebLicenseInstaller.loadOELicenseFromWeb();
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        oemolistream ifs = new oemolistream();
        ifs.open("/Users/jfeng1/tmp.sdf");
        OEGraphMol mol = new OEGraphMol();
        Vector<PropertyMolecule> molList = new Vector<PropertyMolecule>();
        while(oechem.OEReadMolecule(ifs,mol)){
            PropertyMolecule p = new PropertyMolecule(mol);
            molList.add(p);
        }

        JPopupMenu.setDefaultLightWeightPopupEnabled(false);
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(new Dimension(1024,768));
        MolViewer3D mviewer = new MolViewer3D();
        frame.getContentPane().add(mviewer);
        for(PropertyMolecule p:molList) {
            mviewer.addLigand(p);
        }
        frame.setVisible(true);

    }

}
