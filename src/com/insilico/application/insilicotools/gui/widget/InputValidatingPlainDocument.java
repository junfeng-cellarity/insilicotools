package com.insilico.application.insilicotools.gui.widget;
import java.awt.Toolkit;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.PlainDocument;

public abstract class InputValidatingPlainDocument extends PlainDocument {

    protected String oldText = "none";

    protected String potentialOldText = "none";
    protected boolean isValidating = false;

    public String getOldText() {
        return oldText.equals("none") ? "" : oldText;
    }

    public void setOldText(String oldText) {
        this.oldText = oldText;
    }

    public String getPotentialOldText() {
        return potentialOldText;
    }

    public void setPotentialOldText(String potentialOldText) {
        this.potentialOldText = potentialOldText;
    }

    public void insertString(int offset, String s, AttributeSet attributeSet) throws BadLocationException {
        StringBuffer sb = null;
        if (!isValidating()) {
            sb = new StringBuffer(s == null ? "" : s);
            char c;
            for (int i = 0; i < sb.length(); i++) {
                c = sb.charAt(i);
                if (isCharacterValid(c, offset)) {
                    potentialOldText = this.getText(0, this.getLength());
                } else {
                    sb.deleteCharAt(i);
                    Toolkit.getDefaultToolkit().beep();
                }
            }

            s = sb.toString();
        }

        super.insertString(offset, s, attributeSet);
    }


    public boolean isValidating() {
        return isValidating;
    }

    public void setValidating(boolean validating) {
        isValidating = validating;
    }

    public abstract boolean isCharacterValid(char c, int offset);
}
