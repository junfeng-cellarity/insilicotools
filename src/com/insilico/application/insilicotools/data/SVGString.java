package com.insilico.application.insilicotools.data;

import com.google.common.base.Strings;

/**
 * Created by jfeng1 on 9/10/15.
 */
public class SVGString {
    String content;

    public SVGString(String svgString) {
        content = svgString;
    }

    public String getContent() {
        return content;
    }

    public boolean isEmpty(){
        return Strings.isNullOrEmpty(content);
    }
}
