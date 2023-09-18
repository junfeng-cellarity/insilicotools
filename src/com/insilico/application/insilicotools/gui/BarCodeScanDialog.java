package com.insilico.application.insilicotools.gui;

import com.insilico.application.insilicotools.gui.util.ImageUtil;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.jidesoft.dialog.ButtonPanel;
import com.jidesoft.dialog.StandardDialog;
import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.util.TempFile;
import org.apache.xmlrpc.XmlRpcException;

import javax.swing.*;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import java.awt.*;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.MalformedURLException;

/**
 * Created by jfeng1 on 1/27/16.
 */
public class BarCodeScanDialog extends StandardDialog{
    private JTextArea textArea;
    JButton okBtn;
    JButton copyBtn;
    JButton excelBtn;
    private Timer timer;
    private static boolean updating = false;

    public BarCodeScanDialog() {
        super();
        timer = new Timer(1000, new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                checkBarCode();
            }
        });

        textArea = new JTextArea();
        textArea.getDocument().addDocumentListener(new DocumentListener() {

            public void insertUpdate(DocumentEvent e) {
                if(timer.isRunning()){
                    timer.stop();
                }
                timer.start();
            }

            @Override
            public void removeUpdate(DocumentEvent e) {

            }

            @Override
            public void changedUpdate(DocumentEvent e) {

            }
        });
        okBtn = new JButton("Close");
        okBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                setVisible(false);
            }
        });

        copyBtn = new JButton("Copy Text");
        copyBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(textArea.getText().isEmpty()){
                    JOptionPane.showMessageDialog(BarCodeScanDialog.this,"No barcode available.");
                    return;
                }
                getToolkit().getSystemClipboard().setContents(new StringSelection(textArea.getText()),null);
            }
        });

        excelBtn = new JButton("Open in excel");
        excelBtn.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(textArea.getText().isEmpty()){
                    JOptionPane.showMessageDialog(BarCodeScanDialog.this,"No barcode available.");
                    return;
                }
                String content = textArea.getText();
                String[] lines = content.split("\n");
                int rowNum = lines.length;
                try {
                    File tempFile = TempFile.createTempFile("barcode", ".xls");
                    HSSFWorkbook excelWorkBook = new HSSFWorkbook();
                    HSSFSheet sheet1 = excelWorkBook.createSheet();
                    for (int i = 0; i < rowNum; i++) {
                        String line = lines[i];
                        String[] args = line.split("\\s+");
                        int colNum = Math.max(2,args.length);
                        if(colNum==0) {
                            continue;
                        }

                        HSSFRow row = sheet1.createRow(i+1);
                        for (int j = 0; j < colNum; j++) {
                            HSSFCell cell = row.createCell(j+1);
                            cell.setCellValue(args[j]);
                        }
                    }
                    excelWorkBook.write(new FileOutputStream(tempFile));
                    System.out.println(tempFile.getAbsolutePath());
                    Desktop.getDesktop().open(tempFile);
                } catch (IOException e1) {
                    e1.printStackTrace();
                    JOptionPane.showMessageDialog(BarCodeScanDialog.this,e1.getMessage());
                }
            }
        });

    }

    private void checkBarCode(){
        if(updating){
            return;
        }
        try {
            final StringBuilder sb = new StringBuilder();
            String text = textArea.getText();
            if(text.trim().isEmpty()){
                return;
            }
            String[] lines = text.split("\n");
            for(String line:lines){
                String[] args = line.split("\\s+");
                if(args.length==1){
                    String barcode = args[0].trim();
                    String result = ChemFunc.getTareWeight(barcode);
                    sb.append(String.format("%s %s\n",barcode,result));
                }else{
                    sb.append(String.format("%s\n",line));
                }
            }
            updating = true;
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    textArea.setText(sb.toString());
                }
            });
            updating = false;

        } catch (MalformedURLException e) {
            e.printStackTrace();
        } catch (XmlRpcException e) {
            e.printStackTrace();
        }
    }

    @Override
    public JComponent createBannerPanel() {
        ImageIcon imageIcon = ImageUtil.resizeIcon(new ImageIcon(getClass().getClassLoader().getResource("Barcode-Scanner-icon.png")));
        JLabel label = new JLabel("Scan Your Vial One By One", imageIcon, SwingConstants.CENTER);
        return label;
    }

    @Override
    public JComponent createContentPanel() {
        return textArea;
    }

    @Override
    public ButtonPanel createButtonPanel() {
        ButtonPanel p = new ButtonPanel();
        p.setAlignment(SwingConstants.CENTER);
        p.addButton(copyBtn);
        p.addButton(okBtn);
        p.addButton(excelBtn);
        return p;
    }

    public static void main(String[] args) {
        BarCodeScanDialog dialog = new BarCodeScanDialog();
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.setSize(new Dimension(800,600));
        dialog.setVisible(true);
    }

}
