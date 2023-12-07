package com.insilico.application.insilicotools;

import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import ch.qos.logback.classic.LoggerContext;
import chemaxon.license.LicenseHandler;
import chemaxon.license.LicenseManager;
import chemaxon.license.LicenseProcessingException;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.insilico.application.insilicotools.util.OEChemWebLicenseInstaller;
import com.insilico.application.insilicotools.util.OSValidator;
import com.jgoodies.looks.plastic.Plastic3DLookAndFeel;
import com.jgoodies.looks.plastic.PlasticLookAndFeel;
import com.jgoodies.looks.plastic.theme.*;
import org.json.simple.JSONObject;
import org.slf4j.LoggerFactory;

import javax.swing.*;
import javax.swing.plaf.metal.MetalLookAndFeel;
import java.awt.*;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import com.google.common.io.Files;
import openeye.oechem.*;
import java.net.URL;
import java.util.TimeZone;

/**
 * Created by jfeng1 on 9/4/15.
 */
public class inSilicoTools extends JFrame{
//    static { System.setProperty("logback.configurationFile", "http://10.74.2.128:8080/insilico/logback.xml");}
    static inSilicoTools _this;
    public inSilicoTools() throws HeadlessException {
        super(new String("InSilico Tools"));
    }

    public static inSilicoTools getInstance(){
        if(_this==null){
            _this = new inSilicoTools();
        }
        return _this;
    }
//    public Logger getLogger(){
//        Logger logger = LoggerFactory.getLogger(inSilicoTools.class);
//        return logger;
//    }

    public void logStartup(){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        String os =System.getProperty("os.name");
        String arch = System.getProperty("os.arch");
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("os",os);
        object.put("arch",arch);
        object.put("function","inSilicoTools");
//        logger.info(object.toJSONString());

    }

    public void logEnumeration(){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("function","Enumeration");
//        logger.info(object.toJSONString());
    }

    public void logCompoundTracking(String project){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("function","Compound Tracking");
        object.put("project",project);
//        logger.info(object.toJSONString());
    }

    public void logDocking(String structureName, int numMols){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("function","Docking");
        object.put("StructureName",structureName);
        object.put("nMols",numMols);
//        logger.info(object.toJSONString());
    }

    public void logSuperposition(String templateSmiles, int numMols){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("function","Superposition");
        object.put("template",templateSmiles);
        object.put("nMols",numMols);
//        logger.info(object.toJSONString());
    }

    public void logPrediction(String properties,int numMols){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("function","Model Prediction");
        object.put("properties",properties);
        object.put("nMols",numMols);
//        logger.info(object.toJSONString());
    }

    public void logConformationGeneration(String smiles){
//        Logger logger = getInstance().getLogger();
        String user = InSlilicoPanel.getInstance().getUserName();
        JSONObject object = new JSONObject();
        object.put("user",user);
        object.put("molecule",smiles);
        object.put("function","ConfGen");
//        logger.info(object.toJSONString());
    }


    public static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                TimeZone.setDefault(TimeZone.getTimeZone("EST5EDT"));
                LoggerContext loggerContext = (LoggerContext) LoggerFactory.getILoggerFactory();
                Logger rootLogger = loggerContext.getLogger("chemaxon");
                rootLogger.setLevel(Level.OFF);

                JPopupMenu.setDefaultLightWeightPopupEnabled(false);
                ToolTipManager.sharedInstance().setLightWeightPopupEnabled(false);
                if(System.getProperty("os.name").equals("Linux")){
                    try {
                        UIManager.setLookAndFeel(new MetalLookAndFeel());
                    } catch (UnsupportedLookAndFeelException e) {
                        e.printStackTrace();
                    }
                }else {
                    try {
                        PlasticLookAndFeel.setPlasticTheme(new SkyBlue());
                        UIManager.setLookAndFeel(new Plastic3DLookAndFeel());
                    } catch (UnsupportedLookAndFeelException e) {
                        e.printStackTrace();
                    }
                }

                try {
//                    OEChemWebLicenseInstaller.loadOELicenseFromWeb();
                    oechem.OEAddLicenseFromHttp(new URL("http://10.74.2.128:8080/insilico_tools/oe_license.txt"));
//                    LicenseManager.setLicenseFile("/home/jfeng/.chemaxon/license.cxl");
                    LicenseManager.setLicenseFile("http://10.74.2.128:8080/insilico_tools/license.cxl");
//                    LicenseHandler.getInstance().checkLicense("Marvin Applets");
                    LicenseHandler.getInstance().checkLicense("Marvin Beans");
//                    LicenseHandler.getInstance().checkLicense("JChem Base");
//                    LicenseHandler.getInstance().checkLicense("Protonation Plugin Group");
                } catch (IOException e) {
                    JOptionPane.showMessageDialog(null,e.getMessage());
                    return;
                } catch (LicenseProcessingException e) {
                    JOptionPane.showMessageDialog(null,e.getMessage());
                    return;
                }

                if(OSValidator.isUnix()) {
                    inSilicoTools.getInstance().loadLibraryFromResource("libGraphMolWrap.so");
                    inSilicoTools.getInstance().loadLibraryFromResource("libGLEW.so.1.5");
                    inSilicoTools.getInstance().loadLibraryFromResource("libgluegen-rt.so");
                    inSilicoTools.getInstance().loadLibraryFromResource("libjogl_desktop.so");
//                    inSilicoTools.getInstance().loadLibraryFromResource("libnativewindow_awt.so");
                    inSilicoTools.getInstance().loadLibraryFromResource("libnativewindow_x11.so");
                    inSilicoTools.getInstance().loadLibraryFromResource("libnewt.so");
                    inSilicoTools.getInstance().loadLibraryFromResource("libjymol.so");
                }else if(OSValidator.isMac()){
                    inSilicoTools.getInstance().loadLibraryFromResource("libGraphMolWrap.jnilib");
                }else{
                    inSilicoTools.getInstance().loadLibraryFromResource("boost_serialization-vc140-mt-1_65_1.dll");
                    inSilicoTools.getInstance().loadLibraryFromResource("GraphMolWrap.dll");
                }

                JFrame frame = inSilicoTools.getInstance();
//                inSilicoTools.getInstance().logStartup();
                InSlilicoPanel p = InSlilicoPanel.getInstance();
                frame.getContentPane().add(p);
                frame.setJMenuBar(p.getMenuBar());
                frame.setSize(new Dimension(1280,800));
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                frame.setLocationRelativeTo(null);
                frame.setVisible(true);
            }
        });
    }

    private void loadLibraryFromResource(String libName){
        try {
            String[] args = libName.split("\\.");
            InputStream in = getClass().getClassLoader().getResourceAsStream(libName);
            byte[] buffer = new byte[1024];
            int read = -1;
            System.out.println(args[0]);
            File temp = new File(Files.createTempDir(),libName);
            //File temp = File.createTempFile(args[0], "."+args[1]);
            FileOutputStream fos = new FileOutputStream(temp);
            while((read = in.read(buffer)) != -1) {
                fos.write(buffer, 0, read);
            }
            fos.close();
            in.close();

            System.load(temp.getAbsolutePath());
        } catch (IOException e) {
            JOptionPane.showMessageDialog(InSlilicoPanel.getInstance(),e.getMessage());
            e.printStackTrace();
        }
    }
}
