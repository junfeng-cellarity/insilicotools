package com.insilico.application.insilicotools.gui.util;

import javax.swing.*;
import java.awt.*;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URLEncoder;

/**
 * Created by jfeng1 on 10/27/15.
 */
public class ImageUtil {
    private static final String SMI_SERVER = "javelin.corp.My.com";
    private static final int SMI_SERVER_PORT = 5000;


    public static String smi2ImgUrl(String smiles){
        if(smiles==null){
            return "<html><body><p>No molecule available.</p></body></html>";
        }
        String[] args = smiles.split(" ",2);
        String name = null;
        String smi = args[0];
        if(args.length==2){
            name = args[1];
        }
        try {
            URI uri = new URI("http", null, SMI_SERVER, SMI_SERVER_PORT, "/smiles/"+ URLEncoder.encode(smi,"UTF-8"),null,null);
            if(name!=null) {
                return String.format("<html><<body><img src=\"%s\" width=\"300\" height=\"300\">\n<div align=\"center\">%s</div></body></html>", uri.toString(), name);
            }else{
                return String.format("<html><<body><img src=\"%s\" width=\"300\" height=\"300\">\n</body></html>", uri.toString());            }
        } catch (URISyntaxException e) {
            e.printStackTrace();
        } catch (UnsupportedEncodingException e) {
            e.printStackTrace();
        }
        return "";
    }

    public static ImageIcon resizeIcon(ImageIcon icon){
        Image img = icon.getImage() ;
        Image newimg = img.getScaledInstance( 24, 24,  java.awt.Image.SCALE_SMOOTH ) ;
        return new ImageIcon( newimg );
    }

}
