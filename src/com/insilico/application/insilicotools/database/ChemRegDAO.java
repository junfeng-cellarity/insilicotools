package com.insilico.application.insilicotools.database;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemWebLicenseInstaller;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;

import java.io.*;
import java.sql.*;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created by jfeng1 on 9/19/16.
 */
public class ChemRegDAO {
    final String dbURL = "jdbc:oracle:thin:@//10.2.129.34:1725/PRELNR";
    final static String user = "chem_user";
    final static String password = "Myidec123";
    Connection connection;
    PreparedStatement preparedStatement;
    ResultSet rs;

    public static final int SINGLE_BATCH = 1;
    public static final int SMALL_BATCH = 4;
    public static final int MEDIUM_BATCH = 11;
    public static final int LARGE_BATCH = 51;

    static ChemRegDAO _this = null;

    public static ChemRegDAO getInstance(){
        if(_this==null){
            _this = new ChemRegDAO();
        }
        return _this;
    }

    private ChemRegDAO() {
        try {
            Class.forName("oracle.jdbc.OracleDriver");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }
}
