package com.insilico.application.insilicotools.data;

import com.google.common.collect.BoundType;
import com.google.common.collect.Range;
import com.google.common.hash.HashCode;

import java.text.DecimalFormat;

/**
 * Created by jfeng1 on 9/11/15.
 */
public class MolProperty implements Comparable<MolProperty>{
    String property;
    double value;
    boolean isNumerical = true;
    double stddev = -1;
    DecimalFormat df = new DecimalFormat("0.##");
    Range<Double> range = null;
    String unit = "";

    public MolProperty(String property) {
        this.property = property;
        try {
            value = Double.parseDouble(property);
        } catch (NumberFormatException e) {
            isNumerical = false;
        }
    }

    public String getUnit() {
        return unit;
    }

    public void setUnit(String unit) {
        this.unit = unit;
    }

    public void setRange(double lower, double upper){
        range = Range.range(lower, BoundType.OPEN,upper,BoundType.OPEN);
    }

    public void setStddev(double stddev) {
        this.stddev = stddev;
    }

    public MolProperty(boolean isNumerical) {
        this.isNumerical = isNumerical;
    }

    public String getProperty() {
        return property;
    }

    public double getValue() throws NumberFormatException{
        if(isNumerical) {
            return value;
        }else{
            throw new NumberFormatException("Not a numerical value.");
        }
    }

    public boolean isNumerical() {
        return isNumerical;
    }

    @Override
    public String toString() {
        if(isNumerical){
            if(stddev>0){
                return String.format("%s \u00B1 %s%s",df.format(value),df.format(stddev),unit);
            }
//            else if(range!=null){
//                return String.format("%s (%5.2f ~ %5.2f)",df.format(value),range.lowerEndpoint(),range.upperEndpoint());
//            }
            else{
                return df.format(value)+unit;
            }
        }else{
            return property;
        }
    }

    @Override
    public int hashCode() {
        return HashCode.fromString(property).hashCode();
    }

    @Override
    public int compareTo(MolProperty o) {
        if(this.isNumerical&&o.isNumerical){
            return new Double(this.getValue()).compareTo(new Double(o.getValue()));
        }else{
            return property.compareTo(o.getProperty());
        }
    }
}
