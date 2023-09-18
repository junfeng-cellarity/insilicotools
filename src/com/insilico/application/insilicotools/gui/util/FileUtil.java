package com.insilico.application.insilicotools.gui.util;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;
import java.io.*;

/**
 * Created by jfeng1 on 9/26/15.
 */
public class FileUtil {

    public static void saveToFile(File currentDir, FileNameExtensionFilter fileNameExtensionFilter, FileFunctor functor){
        JFileChooser fc = new JFileChooser(currentDir);
        fc.addChoosableFileFilter(fileNameExtensionFilter);
        int option = fc.showSaveDialog(null);
        if(option!=JFileChooser.APPROVE_OPTION){
            return;
        }

        File selectedFile = fc.getSelectedFile();
        if(!selectedFile.getAbsolutePath().endsWith(fileNameExtensionFilter.getExtensions()[0])){
            selectedFile = new File(String.format("%s.%s", selectedFile.getAbsolutePath(), fileNameExtensionFilter.getExtensions()[0]));
        }
        if (selectedFile.exists()) {
            option = JOptionPane.showConfirmDialog(null, "File exists, overwrite?", "Confirmation", JOptionPane.YES_NO_OPTION);
            if(option== JOptionPane.YES_OPTION){
                functor.execute(selectedFile);
            }
        } else {
            functor.execute(selectedFile);
        }
    }

    public static String readFileToString(File file) throws IOException {
        if(file!=null&&file.exists()&&file.canRead()){
            BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
            StringBuilder sb = new StringBuilder();
            String s;
            while((s = reader.readLine())!=null){
                sb.append(s);
                sb.append("\n");
            }
            reader.close();
            return sb.toString();
        }else{
            throw new IOException("File not available.");
        }
    }

    public static void writeStringToFile(File file, String string) throws IOException{
        if((file.exists()&&file.canWrite())||!file.exists()){
            BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file)));
            writer.write(string);
            writer.close();
        }else{
            throw new IOException("File not writable.");
        }
    }
}
