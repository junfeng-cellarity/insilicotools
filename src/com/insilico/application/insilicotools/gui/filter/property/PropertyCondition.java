package com.insilico.application.insilicotools.gui.filter.property;

public class PropertyCondition {
    public static final String GREATER_THAN = ">";
    public static final String LESS_THAN = "<";
    public static final String EQUALS = "=";
    public static final String BETWEEN = "between";

    private final String property;
    private final String operator;
    private final String text;
    private final String text2;


    public PropertyCondition(String property, String operator, String text,  String text2) {
        this.property = property;
        this.operator = operator;
        this.text = text;
        this.text2 = text2;
    }


    public String getProperty() {
        return property;
    }

    public String getOperator() {
        return operator;
    }

    public String getText() {
        return text;
    }

    public String getText2() {
        return text2;
    }
}