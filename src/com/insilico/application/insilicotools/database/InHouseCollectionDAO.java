package com.insilico.application.insilicotools.database;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;

import java.sql.*;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created by jfeng1 on 9/22/15.
 */
public class InHouseCollectionDAO {
    final String dbURL = "jdbc:postgresql://10.74.2.128/cellarity";
    final String dbURL_virtual = "jdbc:postgresql://10.74.2.128/cellarity_virtual";
    final static String user = "medchem";
    final static String password = "medchem";
    Connection connection;
    PreparedStatement preparedStatement;
    ResultSet rs;

    static InHouseCollectionDAO _this = null;
    public static final int SINGLE_BATCH = 1;
    public static final int SMALL_BATCH = 4;
    public static final int MEDIUM_BATCH = 11;
    public static final int LARGE_BATCH = 51;



    public static InHouseCollectionDAO getInstance(){
        if(_this==null){
           _this = new InHouseCollectionDAO();
        }
        return _this;
    }

    private InHouseCollectionDAO() {
        try {
            Class.forName("org.postgresql.Driver");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void cleanup(Connection connection, ResultSet rs, PreparedStatement preparedStatement) {
        try {
            if(rs!=null) {
                rs.close();
            }
        }catch (SQLException ex) {
            rs = null;
        }
        try {
            if(preparedStatement!=null) {
                preparedStatement.close();
            }
        } catch (SQLException ex) {
            preparedStatement = null;
        }
        try {
            if(connection!=null) {
                connection.close();
            }
        } catch (SQLException ex) {
            connection = null;
        }
    }


    public synchronized int getSimilarMoleculeCountFromCellarity(String smiles, float similarity, int maxHits) throws SQLException{
        if(maxHits>0){
            return maxHits;
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        int count = 0;
        try {
            connection = DriverManager.getConnection(dbURL, user, password);
            String query = "select count(*) from chemicals " +
                    "inner join rdk.fps on chemicals.chemical_id = fps.chemical_id " +
                    "and tanimoto_sml(mfp2,morganbv_fp(mol_from_smiles(?::cstring)))>?;";

            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smiles);
            preparedStatement.setFloat(2,similarity);
            rs = preparedStatement.executeQuery();
            if (rs.next()) {
                count = rs.getInt(1);
            }
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
        return count;
    }


    public synchronized Vector<PropertyMolecule> getSimilarMoleculesFromCellarity(String smiles, float similarity, int maxHits) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
        try {
            connection = DriverManager.getConnection(dbURL, user, password);
            String query = "select chemicals.molfile,chemicals.corporate_id from chemicals " +
                    "inner join rdk.fps on chemicals.chemical_id = fps.chemical_id " +
                    "and tanimoto_sml(mfp2,morganbv_fp(mol_from_smiles(?::cstring)))>?";
            if(maxHits>0){
                query = query+" limit ?;";
            }

            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smiles);
            preparedStatement.setFloat(2,similarity);
            if(maxHits>0) {
                preparedStatement.setInt(3, maxHits);
            }
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                String corporate_id = rs.getString("corporate_id");
                OEGraphMol mol = ChemFunc.getMolFromMolString(molfile,OEFormat.SDF);
                mol.SetTitle(corporate_id);
                oechem.OESetSDData(mol,"CY_Number",corporate_id);
                propertyMolecules.add(new PropertyMolecule(mol));
            }
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
        return propertyMolecules;
    }


    public synchronized int getCountFromSmartsFromCellarity(String smarts, int limit) throws SQLException {
        try {
            connection = DriverManager.getConnection(dbURL, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(chemical_id) from rdk.mols where m@>?::qmol limit %d", limit);
            }else{
                query = "select count(chemical_id) from rdk.mols b where m@>?::qmol";
            }
            System.out.println(query);
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            int count = 0;
            if (rs.next()) {
                count = rs.getInt(1);
            }
            return count;
        } finally {
            try {
                if(rs!=null) {
                    rs.close();
                }
            }catch (java.sql.SQLException ex) {
                rs = null;
            }
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
    }


    public synchronized Vector<String> getMolFromSmartsFromCellarity(String smarts, int limit) throws SQLException {
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL, user, password);
            String query;
            if(limit>0) {
                query = String.format("select a.orig_molfile,a.corporate_id from chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id "+
                        " and m@>?::qmol limit %d", limit);
            }else{
                query = "select a.orig_molfile,a.corporate_id from chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        " and m@>?::qmol";
            }
            System.out.println(query);
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("orig_molfile");
                String corporate_id = rs.getString("corporate_id");
                OEGraphMol molFromMolString = ChemFunc.getMolFromMolString(molfile, OEFormat.SDF);
                molFromMolString.SetTitle(corporate_id);
                molFiles.add(ChemFunc.getMolString(molFromMolString));
            }
            return molFiles;
        }catch(Exception e){
            e.printStackTrace();
            return null;
        } finally {
            try {
                if(rs!=null) {
                    rs.close();
                }
            }catch (java.sql.SQLException ex) {
                rs = null;
            }
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
    }

    public synchronized void setChiralSSS(boolean chiralSSS, String db_url) throws SQLException{
        try{
            connection = DriverManager.getConnection(db_url, user, password);
            if(chiralSSS) {
                preparedStatement = connection.prepareStatement("set rdkit.do_chiral_sss=true");
            }else{
                preparedStatement = connection.prepareStatement("set rdkit.do_chiral_sss=false");
            }
            preparedStatement.executeUpdate();
        }finally {
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
    }

    public synchronized String getCYNumberFromMol(String molfile) throws SQLException{
        setChiralSSS(true, dbURL);
        try {
            String cy_number = null;
            connection = DriverManager.getConnection(dbURL, user, password);
            preparedStatement = connection.prepareStatement("select corporate_id from chemicals a inner join rdk.mols b on a.chemical_id=b.chemical_id where m@=mol_from_ctab(?::cstring)");
            preparedStatement.setString(1,molfile);
            rs = preparedStatement.executeQuery();
            if (rs.next()) {
                cy_number = rs.getString("corporate_id");
            }
            return cy_number;
        } finally {
            try {
                rs.close();
            }catch (java.sql.SQLException ex) {
                rs = null;
            }
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
    }

    public synchronized String getCYTautomerFromMol(String molFile) throws SQLException{
        setChiralSSS(false,dbURL);
        try {
            String cy_number = null;
            connection = DriverManager.getConnection(dbURL, user, password);
            preparedStatement = connection.prepareStatement("select corporate_id from chemicals a inner join rdk.mols b on a.chemical_id=b.chemical_id where m@>mol_from_ctab(?::cstring) and m<@mol_from_ctab(?::cstring)");
            preparedStatement.setString(1,molFile);
            preparedStatement.setString(2,molFile);
            rs = preparedStatement.executeQuery();
            if (rs.next()) {
                cy_number = rs.getString("corporate_id");
            }
            return cy_number;
        } finally {
            try {
                rs.close();
            }catch (java.sql.SQLException ex) {
                rs = null;
            }
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
    }

    public synchronized Vector<PropertyMolecule> getMolFromCYNumberBatch(Vector<String> cyNumberList, ProgressReporter progressReporter) throws SQLException {
        Vector<PreparedStatement> statements = new Vector<PreparedStatement>();
        Vector<ResultSet> resultSets = new Vector<ResultSet>();
        Vector<PropertyMolecule> molecules = new Vector<PropertyMolecule>();
        HashMap<String,PropertyMolecule> dict = new HashMap<String, PropertyMolecule>();
        Vector<String> bioList = new Vector<>();
        bioList.addAll(cyNumberList);
        try {
            connection = DriverManager.getConnection(dbURL, user, password);
            int totalNumberOfValuesLeftToBatch = cyNumberList.size();
            final int total =cyNumberList.size();
            while(totalNumberOfValuesLeftToBatch>0){
                int batchSize = SINGLE_BATCH;
                if ( totalNumberOfValuesLeftToBatch >= LARGE_BATCH ) {
                    batchSize = LARGE_BATCH;
                } else if ( totalNumberOfValuesLeftToBatch >= MEDIUM_BATCH ) {
                    batchSize = MEDIUM_BATCH;
                } else if ( totalNumberOfValuesLeftToBatch >= SMALL_BATCH ) {
                    batchSize = SMALL_BATCH;
                }
                totalNumberOfValuesLeftToBatch -= batchSize;
                if(progressReporter!=null){
                    progressReporter.reportProgress("Loading molecules from CY Numbers...", 100*(total-totalNumberOfValuesLeftToBatch)/total);
                }
                StringBuilder inClause = new StringBuilder();
                boolean firstValue = true;
                for (int i=0; i < batchSize; i++) {
                    if ( firstValue ) {
                        inClause.append('?');
                        firstValue = false;
                    } else {
                        inClause.append(',');
                        inClause.append('?');
                    }
                }
                PreparedStatement stmt = connection.prepareStatement("select orig_molfile, corporate_id, chemical_id from chemicals " +
                        "where corporate_id in (" + inClause.toString() + ')');

                for(int i=1;i<=batchSize;i++){
//                    System.out.println(cyNumberList.get(i-1));
                    String cy_number = cyNumberList.get(i - 1);
                    stmt.setString(i, cy_number);
                }

                for(int i=batchSize-1;i>=0;i--){
                    String removed = cyNumberList.remove(i);
                }
                statements.add(stmt);
                ResultSet rs = stmt.executeQuery();
                resultSets.add(rs);
                while(rs.next()){
                    String molFile = rs.getString("orig_molfile");
                    String cy_number = rs.getString("corporate_id");
                    int chemical_id = rs.getInt("chemical_id");
                    OEGraphMol mol = ChemFunc.getMolFromMolString(molFile, OEFormat.SDF);
                    mol.SetTitle(cy_number);
                    PropertyMolecule propertyMolecule = new PropertyMolecule(mol);
                    propertyMolecule.addProperty("cdd_id",""+chemical_id);
                    dict.put(cy_number,propertyMolecule);
                }
                rs.close();
                stmt.close();
            }
            for(String cy_number:bioList){
                if(dict.containsKey(cy_number)){
                    molecules.add(dict.get(cy_number));
                }
            }
        } finally {
            for(ResultSet rs:resultSets){
                try {
                    if(rs!=null) {
                        rs.close();
                    }
                } catch (SQLException e) {
                    e.printStackTrace();
                }
            }

            for(PreparedStatement ps:statements){
                try {
                    if(ps!=null) {
                        ps.close();
                    }
                } catch (SQLException e) {
                    e.printStackTrace();
                }
            }

            if(connection!=null){
                try {
                    connection.close();
                } catch (SQLException e) {
                    connection = null;
                }
            }
        }
        return molecules;
    }


    public static void main(String[] args) {
        InHouseCollectionDAO dao = InHouseCollectionDAO.getInstance();
        try {
            Vector<String> mols = dao.getMolFromSmartsFromCellarity("*N(*)C(=O)c1c[nH]nc1C1CCCN(*)C1", -1);
            for(String m:mols) {
                System.out.println(m);
            }
//            String molFile = dao.getMolFromBioNumber("BIO-0052313");
//            OEGraphMol oemol = new OEGraphMol();
//            oemolistream ifs = new oemolistream();
//            ifs.SetFormat(OEFormat.SDF);
//            ifs.openstring(molFile);
//            oechem.OEReadMolecule(ifs, oemol);
//            ifs.close();
//            String smiles = oechem.OEMolToSmiles(oemol);
//            System.out.println(dao.getBioNumberFromMol(smiles));

        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    public synchronized OEGraphMol getMolFromCYNumber(String cy_number) throws SQLException{
        try {
            String molfile = null;
            int chemical_id = -1;
            connection = DriverManager.getConnection(dbURL, user, password);
            preparedStatement = connection.prepareStatement("select orig_molfile,chemical_id from chemicals where corporate_id = ?");
            preparedStatement.setString(1,cy_number);
            rs = preparedStatement.executeQuery();
            if (rs!=null&&rs.next()) {
                molfile = rs.getString("orig_molfile");
                chemical_id = rs.getInt("chemical_id");
            }
            if(molfile!=null) {
                OEGraphMol mol = ChemFunc.getMolFromMolString(molfile, OEFormat.MDL);
                if (mol != null) {
                    oechem.OESetSDData(mol, "cdd_id", "" + chemical_id);
                    mol.SetTitle(cy_number);
                    return mol;
                }
            }
        } finally {
            try {
                if(rs!=null) {
                    rs.close();
                }
            }catch (java.sql.SQLException ex) {
                rs = null;
            }
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
        return null;
    }

    public synchronized String getVIRNumberFromMol(String mol_string) throws SQLException {
        setChiralSSS(true,dbURL_virtual);
        try {
            String cy_number = null;
            connection = DriverManager.getConnection(dbURL_virtual, user, password);
            preparedStatement = connection.prepareStatement("select corporate_id from chemicals a inner join rdk.mols b on a.chemical_id=b.chemical_id where m@=mol_from_ctab(?::cstring)");
            preparedStatement.setString(1,mol_string);
            rs = preparedStatement.executeQuery();
            if (rs.next()) {
                cy_number = rs.getString("corporate_id");
            }
            return cy_number;
        } finally {
            try {
                rs.close();
            }catch (java.sql.SQLException ex) {
                rs = null;
            }
            try {
                preparedStatement.close();
            } catch (java.sql.SQLException ex) {
                preparedStatement = null;
            }
            try {
                connection.close();
            } catch (java.sql.SQLException ex) {
                connection = null;
            }
        }
    }
}
