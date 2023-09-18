package com.insilico.application.insilicotools.gui.util;

import chemaxon.license.LicenseManager;
import chemaxon.license.LicenseProcessingException;
import chemaxon.marvin.beans.MSketchPane;
import chemaxon.marvin.beans.MViewPane;
import chemaxon.marvin.common.UserSettings;
import chemaxon.marvin.paint.DispOptConsts;

import javax.swing.*;
import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Created by jfeng1 on 12/31/15.
 */
public class MarvinFactory {
    private static MSketchPane sketchPane;
    public static MSketchPane getCompoundSketcher(){
        if(sketchPane==null){
            sketchPane = new MSketchPane();
            sketchPane.setSketchMode(MSketchPane.SM_SELECT_LASSO);
            UserSettings settings = new UserSettings();
            settings.setSketchMolbg2d(Color.white);
            settings.setSketchRendering2d("WIREFRAME");
            settings.setProperty("menuConfig","http://javelin.corp.My.com:8080/insilico/chemaxon_style.xml");
            sketchPane.setUserSettings(settings);
        }
        return sketchPane;
    }

    public static MSketchPane createCompoundSketcher(){
        MSketchPane sketchPane = new MSketchPane();
        sketchPane.setSketchMode(MSketchPane.SM_SELECT_LASSO);
        UserSettings settings = new UserSettings();
        settings.setProperty("menuConfig","http://javelin.corp.My.com:8080/insilico/chemaxon_style.xml");
        settings.setSketchMolbg2d(Color.white);
        settings.setSketchRendering2d("WIREFRAME");
        sketchPane.setUserSettings(settings);
        return sketchPane;
    }
    public static MViewPane createViewPane() {
        MViewPane viewPane = new MViewPane();
        viewPane.setMolbg(Color.white);
        viewPane.setRendering(DispOptConsts.WIREFRAME_RENDERING_S);
        return viewPane;
    }

    public static void main(String[] args) {
        try {
            LicenseManager.setLicenseFile("http://javelin.corp.My.com:8080/insilico/license.cxl");
        } catch (LicenseProcessingException e) {
            e.printStackTrace();
            JOptionPane.showMessageDialog(null,e.getMessage());
        }
        JPanel mainPanel = new JPanel(new BorderLayout());
        final MSketchPane sketcher = MarvinFactory.getCompoundSketcher();
        mainPanel.add(sketcher,BorderLayout.CENTER);
        JPanel p = new JPanel();
        JButton btn = new JButton("Paste");
        btn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
                clipboard.setContents(new StringSelection("C#N"),null);
                sketcher.doPaste();
            }
        });
        p.add(btn);
        mainPanel.add(p,BorderLayout.SOUTH);
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.getContentPane().add(mainPanel);
        frame.setSize(new Dimension(800,600));
        frame.setVisible(true);
    }

}
