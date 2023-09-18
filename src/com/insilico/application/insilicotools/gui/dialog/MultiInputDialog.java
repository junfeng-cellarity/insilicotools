package com.insilico.application.insilicotools.gui.dialog;

import com.jidesoft.range.IntegerRange;
import com.jidesoft.range.NumericRange;
import com.jidesoft.range.Range;

import javax.swing.*;
import javax.swing.text.NumberFormatter;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.HashMap;

/**
 * Created by jfeng1 on 8/1/16.
 */
public class MultiInputDialog extends JDialog {
    JButton okBtn;
    JButton cancelBtn;
    HashMap<String,Number> valueMap;
    HashMap<String,Number> resultMap;
    HashMap<String,Range> rangeMap;
    HashMap<String,JFormattedTextField> fieldMap;
    boolean isCommitted = false;

    public MultiInputDialog(JFrame owner, HashMap<String,Number> initialOptions, HashMap<String,Range> rangeMap) {
        super(owner);
        setModal(true);
        valueMap = initialOptions;
        resultMap = new HashMap<String, Number>();
        this.rangeMap = rangeMap;
        fieldMap = new HashMap<String, JFormattedTextField>();
        int size = initialOptions.size();
        JPanel contentPanel = new JPanel(new BorderLayout());
        JPanel listPanel = new JPanel(new GridLayout(size,1));
        for(final String key:initialOptions.keySet()){
            Number value = initialOptions.get(key);
            JPanel p = new JPanel(new GridLayout(1,2));
            JPanel p1 = new JPanel(new FlowLayout(FlowLayout.CENTER));
            Range range = rangeMap.get(key);
            p1.add(new JLabel(String.format("%s (%s,%s):",key, range.minimum(), range.maximum())));
            p.add(p1);
            NumberFormatter nf = new NumberFormatter();
            nf.setMaximum(range.maximum());
            nf.setMinimum(range.minimum());
            final JFormattedTextField field = new JFormattedTextField(nf);
            fieldMap.put(key,field);
            field.setText(value.toString());
            field.setPreferredSize(new Dimension(150,30));
            p.add(field);
            listPanel.add(p);
        }
        listPanel.setBorder(BorderFactory.createTitledBorder("Options"));
        contentPanel.add(new JScrollPane(listPanel),BorderLayout.CENTER);
        JPanel btnPanel = new JPanel();
        okBtn = new JButton("OK");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                for(String key:valueMap.keySet()){
                    JFormattedTextField field = fieldMap.get(key);
                    if(field.getText().isEmpty()){
                        resultMap.put(key,valueMap.get(key));
                    }else{
                        try {
                            resultMap.put(key, NumberFormat.getInstance().parse(field.getText()));
                        } catch (ParseException e1) {
                            resultMap.put(key,valueMap.get(key));
                        }
                    }
                }
                isCommitted = true;
                setVisible(false);
            }
        });
        btnPanel.add(okBtn);

        cancelBtn = new JButton("Cancel");
        cancelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                isCommitted = false;
                setVisible(false);
            }
        });
        btnPanel.add(cancelBtn);
        contentPanel.add(btnPanel,BorderLayout.SOUTH);

        setContentPane(contentPanel);
    }

    public boolean isCommitted() {
        return isCommitted;
    }

    public HashMap<String, Number> getResultMap() {
        return resultMap;
    }

    public static void main(String[] args) {
        HashMap<String,Number> a = new HashMap<String, Number>();
        a.put("test",1);
        a.put("testtest",2.5);
        a.put("testtesttest",3);

        HashMap<String,Range> b = new HashMap<String, Range>();
        b.put("test",new NumericRange(0.5,2));
        b.put("testtest", new IntegerRange(1,3));
        b.put("testtesttest",new IntegerRange(2,4));

        MultiInputDialog dialog = new MultiInputDialog(null, a,b);
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.pack();
        dialog.setVisible(true);

    }
}
