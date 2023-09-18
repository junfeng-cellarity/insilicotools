package com.insilico.application.insilicotools.gui.widget;

import com.insilico.application.insilicotools.gui.util.JyMolUtilities;
import com.schrodinger.jymol.JyMol;


public class MolObject {
    String name;
    String molString;
    String color;
    String style = JyMolUtilities.STICK_STYLE;
    JyMol jymol;
    String format = "mol";
    int hydrogen_option = ModelingClickConsumer.SHOW_NO_HYDROGEN;
    boolean changed = true;
    boolean zoom = false;
    boolean deleted = false;
    boolean showing = true;
    Double glide_score = null;

    public MolObject(String name, String molString, String color, JyMol jymol, String format) {
        this.name = name;
        this.molString = molString;
        this.color = color;
        this.jymol = jymol;
        this.format = format;
    }

    public MolObject(String name, String molString, String color, JyMol jymol) {
        this(name,molString,color,jymol,"mol");
    }

    public MolObject(MolObject oldObj) {
        this.name = oldObj.name;
        this.molString = oldObj.getMolString();
        this.color = oldObj.color;
        this.jymol = oldObj.jymol;
        this.format = oldObj.format;
    }

    public Double getGlide_score() {
        return glide_score;
    }

    public void setGlide_score(Double glide_score) {
        this.glide_score = glide_score;
    }

    public String getMolString() {
        return molString;
    }

    public String getName() {
        return name;
    }


    public void setZoom(boolean zoom) {
        this.zoom = zoom;
        updateDisplay();
    }

    public String getColor() {
        return color;
    }

    public void setMolString(String molString, boolean fireEvent) {
        this.molString = molString;
        changed = true;
        updateDisplay();
        if (jymol != null && fireEvent) {
            jymol.firePropertyChange(name, true, false);
        }
    }

    public void setColor(String color) {
        this.color = color;
        updateDisplay();
    }

    public void setStyle(String style) {
        this.style = style;
        updateDisplay();
    }

    public void setJymol(JyMol jymol) {
        this.jymol = jymol;
    }

    public int getHydrogen_option() {
        return hydrogen_option;
    }

    public void setHydrogen_option(int hydrogen_option) {
        this.hydrogen_option = hydrogen_option;
    }

    public void show(boolean zoom) {

        if(changed||deleted) {
            jymol.cmd.delete(name);
            jymol.cmd.load(molString, "string", format, name, 0, 0, true, true, false, zoom ? -1 : 0); //for jymol 1.7
//            jymol.cmd.load(molString, "string", format, name, 0, -1, true, true, 0, zoom ? -1 : 0); //for jymol 1.8
            changed = false;
            deleted = false;
        }
        showing = true;
        jymol.cmd.hide("everything",name);
        jymol.cmd.show(style, name);
        if(color.equals(JyMolUtilities.COLOR_B_FACTOR)){
            jymol.cmd.spectrum("b","rainbow_rev","(ss h or ss s or ss l) and e. c");
        }else if(color.equals(JyMolUtilities.COLOR_SECONDARY_STRUCTURE)){
            jymol.cmd.color("red","ss h and e. c");
            jymol.cmd.color("yellow","ss s and e. c");
            jymol.cmd.color("green","ss l and e. c");
        }else if(color.equals(JyMolUtilities.COLOR_EP)){

        }else {
            jymol.cmd.color(color, String.format("(%s and e. c)", name));
        }
        switch (hydrogen_option) {
            case ModelingClickConsumer.SHOW_NO_HYDROGEN:
                jymol.cmd.hide(style,"h.");
                break;
            case ModelingClickConsumer.SHOW_POLAR_HYDROGEN:
                jymol.cmd.hide(style, "h. and (elem c extend 1)");
                break;
            case ModelingClickConsumer.SHOW_ALL_HYDROGEN:
                break;
        }
        if(zoom) {
            jymol.cmd.zoom(name);
        }

    }

    public String getHBondName(){
        return String.format("hbond_%s",name);
    }

    public void remove(){
        jymol.cmd.delete(name,false);
        deleted = true;
        showing = false;
    }

    public boolean isShowing() {
        if(!deleted) {
            return showing;
        }
        return false;
    }

    public void undisplay(){
        remove();
//        if(deleted){
//            return;
//        }
//        showing = false;
//        jymol.cmd.hide("everything",name);
    }

    public void updateDisplay() {
        if(showing) {
            show(this.zoom);
        }
    }

    @Override
    public int hashCode() {
        return com.google.common.base.Objects.hashCode(this.name,this.molString);
    }

    @Override
    public boolean equals(Object obj) {
        if(obj==null){
            return false;
        }
        if(getClass()!=obj.getClass()){
            return false;
        }
        final MolObject other = (MolObject)obj;
        return com.google.common.base.Objects.equal(this.name,other.name)&&com.google.common.base.Objects.equal(this.molString,other.molString);
    }
}
