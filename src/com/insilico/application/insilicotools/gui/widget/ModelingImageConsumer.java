package com.insilico.application.insilicotools.gui.widget;

import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.schrodinger.jymol.JyMolImageConsumer;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

/**
 * Created by jfeng1 on 9/23/16.
 */
public class ModelingImageConsumer extends JyMolImageConsumer{
    MolViewer3D molViewer3D;

    public ModelingImageConsumer(MolViewer3D viewer3D) {
        this.molViewer3D = viewer3D;
    }

    @Override
    public void processImage(BufferedImage image) {
        JFileChooser fc = new JFileChooser();
        fc.setFileFilter(new FileNameExtensionFilter("PNG file","png"));
        fc.setCurrentDirectory(InSlilicoPanel.getInstance().getCurrentDirectory());
        int option = fc.showSaveDialog(molViewer3D);
        if(option == JFileChooser.APPROVE_OPTION){
            File selectedFile = fc.getSelectedFile();
            if(selectedFile.exists()){
                int confirm = JOptionPane.showConfirmDialog(molViewer3D, "File exists, overwrite?");
                if(confirm!=JOptionPane.YES_OPTION){
                    return;
                }
            }
            try {
                ImageIO.write(image,"png",selectedFile);
            } catch (IOException e) {
                e.printStackTrace();
                JOptionPane.showMessageDialog(molViewer3D,e.getMessage());
                return;
            }
            JOptionPane.showMessageDialog(molViewer3D,String.format("File saved to %s",selectedFile.getAbsolutePath()));
        }
    }
}
