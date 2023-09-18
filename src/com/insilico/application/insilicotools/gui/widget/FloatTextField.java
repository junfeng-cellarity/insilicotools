package com.insilico.application.insilicotools.gui.widget;
import javax.swing.text.BadLocationException;

public class FloatTextField extends InputValidatingTextField {

    public class FloatDocument extends InputValidatingPlainDocument {

        public boolean isCharacterValid(char c, int offset) {
            boolean isValid = false;

            try {
                isValid = Character.isDigit(c) ||
                        ((c == '-') && (offset == 0) && !getText(0, 1).equals("-")) ||
                        ((c == '.') && (getText(0, getLength()).indexOf('.') == -1));
            } catch (BadLocationException e) {
                e.printStackTrace();
            }

            return isValid;
        }
    }

    public FloatTextField() {
    }

    public FloatTextField(String text) {
        super(text);
    }

    public FloatTextField(int columns) {
        super(columns);
    }

    public FloatTextField(String text, int columns) {
        super(text, columns);
    }

    public InputValidatingPlainDocument createInputValidatingPlainDocument() {
        return new FloatDocument();
    }
}
