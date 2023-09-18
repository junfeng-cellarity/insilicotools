package com.insilico.application.insilicotools.gui.dialog;

import com.insilico.application.insilicotools.gui.widget.ModelingClickConsumer;
import com.insilico.application.insilicotools.gui.widget.MorphingElement;
import com.insilico.application.insilicotools.gui.widget.MorphingFragment;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.*;

/**
 * Created by jfeng1 on 9/14/16.
 */
public class MolBuildingDialog extends JDialog{
    //    public final static int BUILDING_SUBMODE_VIEW = 6;
    //    public final static int BUILDING_SUBMODE_ATTACH = 7;
    //    public final static int BUILDING_SUBMODE_MODIFY = 8;
    //    public final static int BUILDING_SUBMODE_LINK = 9;
    //    public final static int BUILDING_SUBMODE_DELETE = 10;
    //    public final static int BUILDING_SUBMODE_MODIFY_BOND = 11;
    //    public final static int BUILDING_SUBMODE_DELETE_BOND = 12;

    private final static int VIEW = 0;
    private final static int ATTACH = 1;
    private final static int MODIFY = 2;
    private final static int LINK = 3;
    private final static int DELETE = 4;
    private final static int MODIFY_BOND = 5;
    private final static int DELETE_BOND = 6;
    private final static int ROTATE_TORSION = 7;

    String[] functions = {"View","Attach","Modify","Link","Delete","Modify Bond","Delete Bond","Torsion Rotate"};
    String[] helpText = {"",
            "Click on the atom to be attached",
            "Click on the atom to be changed",
            "Click on two atoms to be linked",
            "Click on the atom to be deleted",
            "Click on two atoms to change bond order",
            "Click on two atoms to delete the bond",
            "Click on two atoms to rotate the bond"
    };
    JButton undoBtn = new JButton("Undo");
    JButton redoBtn = new JButton("Redo");
    JTextField infoField = new JTextField(30);
    ModelingClickConsumer clickConsumer;

    public MolBuildingDialog(ModelingClickConsumer consumer) {
        super();
        clickConsumer = consumer;
        JPanel contentPanel = new JPanel(new BorderLayout());
        contentPanel.setBorder(new LineBorder(Color.GREEN));
        JPanel p = new JPanel(new GridLayout(4,2));
        ButtonGroup bg = new ButtonGroup();
        for(int id=0;id<functions.length;id++){
            final ActionPanel actionPanel = new ActionPanel(functions[id],id==0, id);
            final String text= helpText[id];
            final int currentId = id;
            actionPanel.getStatusBtn().addItemListener(new ItemListener() {
                @Override
                public void itemStateChanged(ItemEvent e) {
                    if(actionPanel.getSelected()){
                        infoField.setText(text);
                    }
                }
            });
            bg.add(actionPanel.getStatusBtn());
            p.add(actionPanel);
        }
        contentPanel.add(p,BorderLayout.CENTER);
        JPanel btnPanel = new JPanel(new GridLayout(2,1));
        JPanel p1 = new JPanel();
        p1.add(undoBtn);
        undoBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                clickConsumer.undo();
            }
        });
        p1.add(redoBtn);
        redoBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                clickConsumer.redo();
            }
        });
        btnPanel.add(p1);
        btnPanel.add(infoField);
        infoField.setEditable(false);
        contentPanel.add(btnPanel,BorderLayout.SOUTH);
        setContentPane(contentPanel);
        setModal(false);
        setAlwaysOnTop(true);
        setSize(new Dimension(300,600));
    }

    class ActionPanel extends JPanel implements MouseListener{
        JRadioButton statusBtn = new JRadioButton();
        Color colorSelected = new Color(5,200,5,125);
        Color colorUnselected = new Color(255,255,255);
        int id;
        public ActionPanel(String title, final boolean isSelected, int id) {
            super(new FlowLayout(FlowLayout.CENTER));
            this.id = id;
            TitledBorder border = new TitledBorder(title);
            border.setTitleColor(Color.BLACK);
            setBorder(border);
            setOpaque(true);
            setSelected(isSelected);
            addMouseListener(this);
//            private final static int VIEW = 0;
//            private final static int ATTACH = 1;
//            private final static int MODIFY = 2;
//            private final static int LINK = 3;
//            private final static int DELETE = 4;
//            private final static int MODIFY_BOND = 5;
//            private final static int DELETE_BOND = 6;
//            private final static int ROTATE_TORSION = 7;

            switch(id){
                case VIEW:
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_VIEW);
                                setBackground(colorSelected);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;
                case ATTACH:
                    final JComboBox cb1 = new JComboBox(MorphingFragment.getCommonGroups());
                    cb1.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setSelected(true);
                            MorphingFragment fragment = (MorphingFragment)cb1.getSelectedItem();
                            clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_ATTACH);
                            clickConsumer.setMorphingFragment(fragment);
                        }
                    });
                    add(cb1);
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                setBackground(colorSelected);
                                MorphingFragment fragment = (MorphingFragment)cb1.getSelectedItem();
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_ATTACH);
                                clickConsumer.setMorphingFragment(fragment);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;
                case MODIFY:
                    final JComboBox cb2 = new JComboBox(MorphingElement.getCommonElements());
                    cb2.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setSelected(true);
                            MorphingElement element = (MorphingElement) cb2.getSelectedItem();
                            clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_MODIFY);
                            clickConsumer.setMorphingElement(element);
                        }
                    });
                    add(cb2);
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                setBackground(colorSelected);
                                MorphingElement element = (MorphingElement) cb2.getSelectedItem();
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_MODIFY);
                                clickConsumer.setMorphingElement(element);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;
                case LINK:
                    final JComboBox cb3 = new JComboBox(new String[]{"Auto","0","1","2","3","4","5","6","7","8","9","10"});
                    cb3.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setSelected(true);
                            String selection = (String)cb3.getSelectedItem();
                            clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_LINK);
                            int linker_length = -1;
                            if(!selection.equals("Auto")){
                                linker_length = Integer.parseInt(selection);
                            }
                            clickConsumer.setMorphing_linker_length(linker_length);
                        }
                    });
                    add(cb3);
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                setBackground(colorSelected);
                                String selection = (String)cb3.getSelectedItem();
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_LINK);
                                int linker_length = -1;
                                if(!selection.equals("Auto")){
                                    linker_length = Integer.parseInt(selection);
                                }
                                clickConsumer.setMorphing_linker_length(linker_length);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;

                case DELETE:
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                setBackground(colorSelected);
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_DELETE);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;

                case MODIFY_BOND:
                    final JComboBox cb4 = new JComboBox(new String[]{"single","double","triple","none"});
                    cb4.addActionListener(new ActionListener() {
                        @Override
                        public void actionPerformed(ActionEvent e) {
                            setSelected(true);
                            String bond_type = (String)cb4.getSelectedItem();
                            clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_MODIFY_BOND);
                            clickConsumer.setMorphing_bond_type(bond_type);
                        }
                    });
                    add(cb4);
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                String bond_type = (String)cb4.getSelectedItem();
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_MODIFY_BOND);
                                clickConsumer.setMorphing_bond_type(bond_type);
                                setBackground(colorSelected);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;

                case DELETE_BOND:
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_DELETE_BOND);
                                setBackground(colorSelected);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;

                case ROTATE_TORSION:
                    final JKnob knob = new JKnob(0,Color.LIGHT_GRAY,Color.BLUE);
                    knob.addMouseMotionListener(new MouseMotionListener() {
                        @Override
                        public void mouseDragged(MouseEvent e) {
                            if(!isSelected) {
                                setSelected(true);
                            }
                            clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_ROTATE_BOND);
                            try {
                                clickConsumer.rotateMolecule(180*knob.getAngle()/Math.PI,true);
                            } catch (Exception e1) {
                                e1.printStackTrace();
                                JOptionPane.showMessageDialog(clickConsumer.getMolViewer3D(),e1.getMessage());
                            }
                        }

                        @Override
                        public void mouseMoved(MouseEvent e) {

                        }
                    });
                    add(knob);
                    statusBtn.addItemListener(new ItemListener() {
                        @Override
                        public void itemStateChanged(ItemEvent e) {
                            if(statusBtn.isSelected()){
                                setBackground(colorSelected);
                                clickConsumer.setBuilding_mode(ModelingClickConsumer.BUILDING_SUBMODE_ROTATE_BOND);
                            }else{
                                setBackground(colorUnselected);
                            }
                        }
                    });
                    break;

                default:
                    break;
            }
        }

        public ActionPanel(String title) {
            this(title,false,0);
        }

        public boolean getSelected(){
            return statusBtn.isSelected();
        }

        public void setSelected(boolean selected){
            statusBtn.setSelected(selected);
            if(selected){
                setBackground(colorSelected);
            }else{
                setBackground(colorUnselected);
            }
        }

        public JRadioButton getStatusBtn() {
            return statusBtn;
        }

        @Override
        public void mouseClicked(MouseEvent e) {
            System.out.println("mouse clicked.");
            setSelected(!getSelected());
        }

        @Override
        public void mousePressed(MouseEvent e) {

        }

        @Override
        public void mouseReleased(MouseEvent e) {

        }

        @Override
        public void mouseEntered(MouseEvent e) {

        }

        @Override
        public void mouseExited(MouseEvent e) {

        }
    }

    public static void main(String[] args) {
        final JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        MolBuildingDialog dialog = new MolBuildingDialog(null);
        dialog.addWindowListener(new WindowListener() {
            @Override
            public void windowOpened(WindowEvent e) {

            }

            @Override
            public void windowClosing(WindowEvent e) {
                f.dispose();
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
        dialog.setResizable(false);
        dialog.setVisible(true);
    }
}
