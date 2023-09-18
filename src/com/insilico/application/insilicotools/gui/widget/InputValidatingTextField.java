package com.insilico.application.insilicotools.gui.widget;

import javax.swing.JTextField;
import javax.swing.text.Document;
public abstract class InputValidatingTextField extends JTextField {

    public InputValidatingTextField() {
        setDocument(createInputValidatingPlainDocument());
    }

    public InputValidatingTextField(String text) {
        super(text);
        setDocument(createInputValidatingPlainDocument());
        setText(text);
        getInputValidatingPlainDocument().setOldText(text);
    }

    public InputValidatingTextField(int columns) {
        super(columns);
        setDocument(createInputValidatingPlainDocument());
    }

    public InputValidatingTextField(String text, int columns) {
        super(text, columns);
        setDocument(createInputValidatingPlainDocument());
        setText(text);
        getInputValidatingPlainDocument().setOldText(text);
    }

    private InputValidatingTextField(Document doc, String text, int columns) {
        super(doc, text, columns);
    }

    public InputValidatingPlainDocument getInputValidatingPlainDocument() {
        return (InputValidatingPlainDocument) getDocument();
    }

    public abstract InputValidatingPlainDocument createInputValidatingPlainDocument();

}
