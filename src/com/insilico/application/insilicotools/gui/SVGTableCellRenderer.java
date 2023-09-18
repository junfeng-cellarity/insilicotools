package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.data.SVGString;
import org.apache.batik.dom.svg.SVGDOMImplementation;
import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderInput;
import org.apache.batik.transcoder.TranscodingHints;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.table.DefaultTableCellRenderer;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import org.apache.batik.transcoder.image.ImageTranscoder;
import org.apache.batik.util.SVGConstants;

public class SVGTableCellRenderer extends DefaultTableCellRenderer{
    int lineWidth = 1;
    public SVGTableCellRenderer() {
    }

    public SVGTableCellRenderer(int lineWidth){
        this.lineWidth = lineWidth;
    }

    public static BufferedImage convertSVGStringToImage(String svgString, Rectangle rect) throws TranscoderException {
        MyTranscoder transcoder = new MyTranscoder();
        byte[] svgImageArray = svgString.getBytes();
        TranscoderInput transcoderInput = new TranscoderInput(new ByteArrayInputStream(svgImageArray));
        TranscodingHints hints = new TranscodingHints();
        hints.put(ImageTranscoder.KEY_WIDTH, new Float(rect.getWidth()*0.85));
        hints.put(ImageTranscoder.KEY_HEIGHT, new Float(rect.getHeight()*0.85));
        hints.put(ImageTranscoder.KEY_DOM_IMPLEMENTATION, SVGDOMImplementation.getDOMImplementation());
        hints.put(ImageTranscoder.KEY_DOCUMENT_ELEMENT_NAMESPACE_URI, SVGConstants.SVG_NAMESPACE_URI);
        hints.put(ImageTranscoder.KEY_DOCUMENT_ELEMENT, SVGConstants.SVG_SVG_TAG);
        hints.put(ImageTranscoder.KEY_XML_PARSER_VALIDATING, new Boolean(false));
        transcoder.setTranscodingHints(hints);
        transcoder.transcode(transcoderInput,null);
        return transcoder.getImage();

    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
        JPanel p = new JPanel(new BorderLayout());
        p.add(this,BorderLayout.CENTER);
        if(value!=null){
            PropertyMolecule molecule = (PropertyMolecule)value;
            JPanel labelPanel = new JPanel();
            labelPanel.add(new JLabel(molecule.getName()));
            p.add(labelPanel,BorderLayout.SOUTH);
            SVGString svgString = molecule.getSvgString();
            if(!svgString.isEmpty()){
                Rectangle rect = table.getCellRect(row, column, true);
                try {
                    BufferedImage molImage = convertSVGStringToImage(svgString.getContent(),rect);
                    ImageIcon imageIcon = new ImageIcon(molImage);
                    this.setIcon(imageIcon);
                    if(molecule.isSelected()){
                        p.setBorder(new LineBorder(Color.RED,lineWidth));
                    }else{
                        p.setBorder(new LineBorder(Color.BLACK,lineWidth));
                    }
                } catch (TranscoderException e) {
                    this.setIcon(new ImageIcon());
                }
            }
        }else{
            this.setIcon(new ImageIcon());
            p.setBorder(new LineBorder(Color.BLACK,0));
        }

        return p;
    }
}
