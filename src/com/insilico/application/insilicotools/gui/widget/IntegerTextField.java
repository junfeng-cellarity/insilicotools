package com.insilico.application.insilicotools.gui.widget;

public class IntegerTextField extends InputValidatingTextField {

    public class IntegerDocument extends InputValidatingPlainDocument {

        public boolean isCharacterValid(char c, int offset) {
            boolean isValid = false;
            if(Character.isDigit(c)) {
                if (offset != 0 || c != '0') {
                    isValid = true;
                }
            }
            return isValid;
        }
    }

    public IntegerTextField() {
    }

    public IntegerTextField(String text) {
        super(text);
    }

    public IntegerTextField(int columns) {
        super(columns);
    }

    public IntegerTextField(String text, int columns) {
        super(text, columns);
    }

    public InputValidatingPlainDocument createInputValidatingPlainDocument() {
        return new IntegerDocument();
    }
}
