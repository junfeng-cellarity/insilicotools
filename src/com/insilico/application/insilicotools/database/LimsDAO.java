package com.insilico.application.insilicotools.database;

import com.insilico.application.insilicotools.data.PropertyMolecule;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.lims.Batch;
import com.insilico.application.insilicotools.gui.lims.LimsMolecule;
import com.insilico.application.insilicotools.util.ChemFunc;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;

import java.sql.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Collections;
import java.util.Vector;

/**
 * Created by jfeng1 on 3/22/17.
 */
public class LimsDAO {
    final String db_url_lims = "jdbc:postgresql://javelin.corp.My.com/lims";

    final static String user = "medchem";
    final static String password = "medchem";

    static LimsDAO _this = null;
    Vector<String> projectNames = new Vector<>();

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

    public static LimsDAO getInstance(){
        if(_this==null){
            _this = new LimsDAO();
        }
        return _this;
    }

    private LimsDAO() {
        try {
            Class.forName("org.postgresql.Driver");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }

    public Vector<LimsMolecule> getMolsNoBatch() throws SQLException{
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        Vector<LimsMolecule> limsMolecules = new Vector<>();
        NumberFormat nf = new DecimalFormat("#.##");
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("select cmpd_id,compound_id,structure,theory_mass,project,scientist,timestamp " +
                    "from lims.compounds where not exists (select lims.batch_to_compounds.cmpd_id from lims.batch_to_compounds where lims.compounds.cmpd_id = lims.batch_to_compounds.cmpd_id)");
            rs = statement.executeQuery();
            while(rs.next()){
                LimsMolecule limsMolecule = getLimsMoleculeFromResultSet(rs, nf);
                limsMolecules.add(limsMolecule);
            }
            return limsMolecules;
        }finally {
            cleanup(connection,rs,statement);
        }
    }


    public Vector<LimsMolecule> getMolsByBatch(int batch_id) throws SQLException{
        if(batch_id==-1){
            return getMolsNoBatch();
        }
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        Vector<LimsMolecule> limsMolecules = new Vector<>();
        NumberFormat nf = new DecimalFormat("#.##");
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("select lims.compounds.cmpd_id,compound_id,structure,theory_mass,project,scientist,lims.compounds.timestamp from lims.compounds inner join lims.batch_to_compounds on lims.compounds.cmpd_id = lims.batch_to_compounds.cmpd_id where lims.batch_to_compounds.batch_id = ?");
            statement.setInt(1,batch_id);
            rs = statement.executeQuery();
            int idx = 1;
            while(rs.next()){
                LimsMolecule limsMolecule = getLimsMoleculeFromResultSet(rs, nf);
                limsMolecules.add(limsMolecule);
            }

            Collections.sort(limsMolecules);
            for(LimsMolecule m:limsMolecules){
                m.generateCompoundId(batch_id,idx++);
            }
            return limsMolecules;
        }finally {
            cleanup(connection,rs,statement);
        }
    }

    LimsMolecule getLimsMoleculeFromResultSet(ResultSet rs, NumberFormat nf) throws SQLException {
        int cmpd_id = rs.getInt("cmpd_id");
        String compound_id = rs.getString("compound_id");
        String structure = rs.getString("structure");
        String theory_mass = rs.getString("theory_mass");
        String project = rs.getString("project");
        String scientist = rs.getString("scientist");
        Long timestamp = rs.getLong("timestamp");
        if(rs.wasNull()){
            timestamp = null;
        }
        OEGraphMol mol = ChemFunc.getMolFromMolString(structure, OEFormat.SDF);
        PropertyMolecule propertyMolecule = new PropertyMolecule(mol);
        LimsMolecule limsMol = new LimsMolecule(propertyMolecule,cmpd_id,compound_id,theory_mass,project,scientist, null,null, timestamp);
        return limsMol;
    }

    public Batch getBatch(int batch_id) throws SQLException{
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("select batch_name, date, chemist from lims.batches where batch_id = ?");
            statement.setInt(1,batch_id);
            rs = statement.executeQuery();
            if(rs.next()){
                String batch_name = rs.getString("batch_name");
                Date date = rs.getDate("date");
                String chemist = rs.getString("chemist");
                return new Batch(batch_id,chemist,batch_name,date);
            }
            return null;
        }finally {
            cleanup(connection,rs,statement);
        }
    }

    public int getNumMolsByBatch(int batch_id){
        if(batch_id==-1){
            try {
                return getMolsNoBatch().size();
            } catch (SQLException e) {
                e.printStackTrace();
            }
            return 0;
        }else{
            try {
                return getMolsByBatch(batch_id).size();
            } catch (SQLException e) {
                e.printStackTrace();
            }
            return 0;
        }
    }

    public int insertBatch(String batch_name, String chemist) throws SQLException{
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("insert into lims.batches (chemist, batch_name) values (?,?)", PreparedStatement.RETURN_GENERATED_KEYS);
            statement.setString(1,chemist);
            statement.setString(2,batch_name);
            int result = statement.executeUpdate();
            if(result == 1){
                ResultSet generatedKeys = statement.getGeneratedKeys();
                if (generatedKeys.next()) {
                    return generatedKeys.getInt(1);
                }
                else {
                    throw new SQLException("Creating batch failed, no ID obtained.");
                }
            }else{
                throw new SQLException("Creating core failed, no ID obtained.");
            }

        }finally {
            cleanup(connection,null,statement);
        }
    }

    public void deleteBatch(int batch_id) throws SQLException{
        if(batch_id<0){
            return;
        }
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("delete from lims.batches where batch_id = ?");
            statement.setInt(1,batch_id);
            statement.executeUpdate();
        }finally {
            cleanup(connection,null,statement);
        }
    }

    public void updateLimsCompound(LimsMolecule pmol) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims, user, password);
            String query;
            query = "update lims.compounds set  compound_id = ?, structure = ?, theory_mass = ?, project = ?, scientist = ? where cmpd_id = ?";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,pmol.getCompound_id());
            preparedStatement.setString(2,pmol.getPropertyMolecule().getSdfStr());
            preparedStatement.setString(3,pmol.getTheory_mass());
            preparedStatement.setString(4,pmol.getProject());
            preparedStatement.setString(5,pmol.getScientist());
            preparedStatement.setInt(6,pmol.getCmpdId());
            preparedStatement.executeUpdate();
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }

    public int insertLimsCompound(LimsMolecule pmol) throws SQLException {
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("insert into lims.compounds " +
                    "(compound_id, structure,theory_mass, project, scientist, timestamp) values (?,?,?,?,?,?)",PreparedStatement.RETURN_GENERATED_KEYS);
            String compound_id = pmol.getCompound_id();
            String structure = ChemFunc.getMolString(pmol.getPropertyMolecule().getMol());
            String theory_mass = pmol.getTheory_mass();
            String scientist = pmol.getScientist();
            String project = pmol.getProject();
            long timestamp = pmol.getTimestamp();
            statement.setString(1, compound_id);
            statement.setString(2, structure);
            statement.setString(3, theory_mass);
            statement.setString(4, project);
            statement.setString(5, scientist);
            statement.setLong(6, timestamp);
            int result = statement.executeUpdate();
            if(result == 1){
                ResultSet generatedKeys = statement.getGeneratedKeys();
                if (generatedKeys.next()) {
                    return generatedKeys.getInt(1);
                }
                else {
                    throw new SQLException("Inserting molecule failed, no ID obtained.");
                }
            }else{
                throw new SQLException("Insert molecule failed, no ID obtained.");
            }
        }finally {
            cleanup(connection,null,statement);
        }
    }



    public void insertCompounds(Vector<LimsMolecule> limsMolecules, ProgressReporter progressReporter){
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            connection.setAutoCommit(false);
            statement = connection.prepareStatement("insert into lims.compounds " +
                    "(compound_id, structure,theory_mass, project, scientist, timestamp) values (?,?,?,?,?,?)");
            int idx = 0;
            for(LimsMolecule pmol:limsMolecules){
                String compound_id = pmol.getCompound_id();
                String structure = ChemFunc.getMolString(pmol.getPropertyMolecule().getMol());
                String theory_mass = pmol.getTheory_mass();
                String project = pmol.getProject();
                String scientist = pmol.getScientist();
                Long timestamp = pmol.getTimestamp();
                statement.setString(1, compound_id);
                statement.setString(2, structure);
                statement.setString(3, theory_mass);
                statement.setString(4, project);
                statement.setString(5,scientist);
                statement.setLong(6,timestamp);
                if(progressReporter!=null){
                    progressReporter.reportProgress("inserting molecule ...",100*idx++/limsMolecules.size());
                }
                statement.addBatch();
            }
            statement.executeBatch();
            connection.commit();

        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,null,statement);
        }
    }

    public void deleteFromBatches(int batch_id, Vector<Integer> cmpdIds) throws SQLException{
        if(batch_id<0){
            return;
        }
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            connection.setAutoCommit(false);
            statement = connection.prepareStatement("delete from lims.batch_to_compounds where batch_id = ? and cmpd_id = ?");
            for(Integer cmpdId:cmpdIds) {
                statement.setInt(1, batch_id);
                statement.setInt(2, cmpdId);
                statement.addBatch();
            }
            statement.executeBatch();
            connection.commit();
        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,null,statement);
        }
    }

    public void deleteCompounds(Vector<Integer> cmpdIds) throws SQLException{
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            connection.setAutoCommit(false);
            statement = connection.prepareStatement("delete from lims.compounds where cmpd_id = ?");
            for(Integer cmpdId:cmpdIds) {
                statement.setInt(1, cmpdId);
                statement.addBatch();
            }
            statement.executeBatch();
            connection.commit();
        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,null,statement);
        }
    }

    public void addCompoundsToBatch(int batch_id, Vector<Integer> cmpdIds) throws SQLException{
        Connection connection = null;
        PreparedStatement statement = null;
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            connection.setAutoCommit(false);
            statement = connection.prepareStatement("insert into lims.batch_to_compounds " +
                    "(batch_id, cmpd_id) values (?,?)");
            for(Integer cmpdId:cmpdIds) {
                statement.setInt(1, batch_id);
                statement.setInt(2, cmpdId);
                statement.addBatch();
            }
            statement.executeBatch();
            connection.commit();
        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,null,statement);
        }
    }


    public static void main(String[] args) {

    }

    public Vector<Batch> getAllBatches() {
        Vector<Batch> batches = new Vector<>();
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        NumberFormat nf = new DecimalFormat("#.##");
        try {
            connection = DriverManager.getConnection(db_url_lims,user,password);
            statement = connection.prepareStatement("select batch_id,date,chemist,batch_name from lims.batches");
            rs = statement.executeQuery();
            while(rs.next()){
                int batch_id = rs.getInt("batch_id");
                Date date = rs.getDate("date");
                String chemist = rs.getString("chemist");
                String batch_name = rs.getString("batch_name");
                batches.add(new Batch(batch_id,chemist,batch_name,date));
            }
            return batches;
        }catch (SQLException e){
            e.printStackTrace();
        }finally {
            cleanup(connection,rs,statement);
        }

        return batches;
    }

    public Vector<String> getMyProjectNames() {
        if(projectNames.isEmpty()) {
            Connection connection = null;
            PreparedStatement statement = null;
            ResultSet rs = null;
            try {
                connection = DriverManager.getConnection(db_url_lims, user, password);
                statement = connection.prepareStatement("select project_name from lims.project order by project_id asc");
                rs = statement.executeQuery();
                while (rs.next()) {
                    String project = rs.getString("project_name");
                    projectNames.add(project);
                }
            } catch (SQLException e) {
                e.printStackTrace();
            } finally {
                cleanup(connection, rs, statement);
            }
        }
        return projectNames;
    }
}
