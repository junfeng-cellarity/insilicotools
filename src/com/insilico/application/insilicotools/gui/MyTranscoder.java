package com.insilico.application.insilicotools.gui;

import org.apache.batik.transcoder.TranscoderException;
import org.apache.batik.transcoder.TranscoderOutput;
import org.apache.batik.transcoder.image.ImageTranscoder;

import java.awt.image.BufferedImage;

/**
 * Created by jfeng1 on 9/8/15.
 */
public class MyTranscoder extends ImageTranscoder{
    BufferedImage image;
    public MyTranscoder() {
        image = null;
    }

    public BufferedImage createImage(int w, int h) {
        image = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
        return image;
    }

    public void writeImage(BufferedImage bufferedImage, TranscoderOutput transcoderOutput) throws TranscoderException {

    }

    public BufferedImage getImage(){
        return image;
    }
}
