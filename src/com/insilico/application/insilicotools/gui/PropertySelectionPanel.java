package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.util.ChemFunc;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Vector;

public class PropertySelectionPanel extends JPanel {
    LinkedHashMap<String,JCheckBox> checkBoxDict;

    public PropertySelectionPanel(String[] properties, String[] adme_models, String[] safety_models, String[] more_models){
        super(new BorderLayout());
        setBorder(new TitledBorder("Properties"));
        JToolBar tbar = new JToolBar();
        JButton selectAllBtn = new JButton("Select all");
        selectAllBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(checkBoxDict!=null){
                    for(JCheckBox b :checkBoxDict.values()){
                        if(Arrays.asList(ChemFunc.TIME_CONSUMING_PROPERTIES_MODELS).contains(b.getText())){
                            continue;
                        }
                        b.setSelected(true);
                    }
                    PropertySelectionPanel.this.firePropertyChange("selection",true,false);
                }
            }
        });
        JButton unSelectAllBtn = new JButton("Unselect all");
        unSelectAllBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (checkBoxDict != null) {
                    for (JCheckBox b : checkBoxDict.values()) {
                        b.setSelected(false);
                    }
                    PropertySelectionPanel.this.firePropertyChange("selection",true,false);
                }
            }
        });

        JButton invertSelectionBtn = new JButton("Invert Selection");
        invertSelectionBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (checkBoxDict != null) {
                    for (JCheckBox b : checkBoxDict.values()) {
                        if (b.isSelected()) {
                            b.setSelected(false);
                        } else {
                            if(Arrays.asList(ChemFunc.TIME_CONSUMING_PROPERTIES_MODELS).contains(b.getText())){
                                continue;
                            }
                            b.setSelected(true);
                        }
                    }
                    PropertySelectionPanel.this.firePropertyChange("selection",true,false);
                }
            }
        });
        tbar.add(selectAllBtn);
        tbar.add(unSelectAllBtn);
        tbar.add(invertSelectionBtn);
        add(tbar, BorderLayout.PAGE_START);

        JPanel p1 = buildPropertyPanel("Properties (RED=Time Consuming)", properties, 3);
        JPanel p2 = buildPropertyPanel("ADME Models", adme_models, 3);
        JPanel p3 = buildPropertyPanel("Safety Models",safety_models,3);
        JPanel p4 = buildPropertyPanel("Miscellaneous",more_models,3);

        JPanel mainPanel = new JPanel();
        BoxLayout boxLayout = new BoxLayout(mainPanel,BoxLayout.Y_AXIS);
        mainPanel.setLayout(boxLayout);
        mainPanel.add(p1);
        mainPanel.add(p2);
        mainPanel.add(p3);
        mainPanel.add(p4);
        add(new JScrollPane(mainPanel), BorderLayout.CENTER);
    }

    private JPanel buildPropertyPanel(String title, String[] properties, int column){
        if(checkBoxDict==null){
            checkBoxDict = new LinkedHashMap<>();
        }
        int nrow = properties.length/2+1;
        JPanel p = new JPanel(new GridLayout(nrow,column));
        p.setBorder(new TitledBorder(title));
        for(String property:properties){
            boolean isDefault = Arrays.asList(ChemFunc.OE_PROPERTIES).contains(property);
            JCheckBox cb = new JCheckBox(property,isDefault);
            if(ChemFunc.PROPERTY_DETAILS.containsKey(property)){
                cb.setToolTipText(ChemFunc.PROPERTY_DETAILS.get(property));
            }
            if(Arrays.asList(ChemFunc.TIME_CONSUMING_PROPERTIES_MODELS).contains(property)){
                cb.setForeground(Color.RED);
            }

            cb.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent e) {
                    PropertySelectionPanel.this.firePropertyChange("selection",true,false);
                    PropertySelectionPanel.this.repaint();
                }
            });
            p.add(cb);
            checkBoxDict.put(property,cb);
        }
        return p;
    }

    public Vector<String> getSelectedProperties(){
        Vector<String> selectedProperties = new Vector<String>();
        for(String p:checkBoxDict.keySet()){
            JCheckBox cb = checkBoxDict.get(p);
            if(cb.isSelected()){
                if(p.equals("MoKa pKa")){
                    selectedProperties.add("MoKa Acidic pKa");
                    selectedProperties.add("MoKa Basic pKa");
                }else if(p.equals("ChemAxon pKa")){
                    selectedProperties.add("ChemAxon Acidic pKa");
                    selectedProperties.add("ChemAxon Basic pKa");
                }else if (p.equals("Human Clearance (route & rate)")){
                    selectedProperties.add("Human Renal Clearance(mL/min/kg)");
                    selectedProperties.add("Human Metabolic Clearance(mL/min/kg)");
                    selectedProperties.add("Clearance Route");
                }
                else {
                    selectedProperties.add(p);
                }
            }
        }
        return selectedProperties;
    }

}
