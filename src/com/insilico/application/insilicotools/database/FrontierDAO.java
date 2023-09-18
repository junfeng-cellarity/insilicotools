package com.insilico.application.insilicotools.database;

import chemaxon.formats.MolExporter;
import chemaxon.struc.Molecule;
import com.insilico.application.insilicotools.data.*;
import com.insilico.application.insilicotools.gui.InSlilicoPanel;
import com.insilico.application.insilicotools.gui.ProgressReporter;
import com.insilico.application.insilicotools.gui.widget.CoreStatus;
import com.insilico.application.insilicotools.util.ChemFunc;
import com.insilico.application.insilicotools.util.OEChemFunc;
import com.google.common.base.Strings;
import openeye.oechem.*;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.sql.*;
import java.util.HashMap;
import java.util.Vector;

/**
 * Created by jfeng1 on 4/20/16.
 */
public class FrontierDAO {
    final String dbURL_Frontier = "jdbc:postgresql://10.74.2.128/frontier";
    final String dbURL_Aldrich = "jdbc:postgresql://10.74.2.128/market_select";
    final String dbURL_GoStar = "jdbc:postgresql://10.74.2.128/GoStar";
    final String dbURL_Enamine = "jdbc:postgresql://10.74.2.128/Enamine";
    final String dbURL_Pharmaron = "jdbc:postgresql://10.74.2.128/Pharmaron";
    final String dbURL_chembl = "jdbc:postgresql://10.74.2.128/chembl_20";
    final String getDbURL_VirtualCmpdLibrary = "jdbc:postgresql://10.74.2.128/cellarity";
    final String getDbURL_CompoundTracker = "jdbc:postgresql://10.74.2.128/compound_tracker";

    final static String user = "medchem";
    final static String password = "medchem";

    static FrontierDAO _this = null;
    HashMap<Integer,String> statusDict = new HashMap<Integer, String>();
    Vector<Status> my_status = new Vector<Status>();
    Vector<Project> my_projects = new Vector<Project>();
    HashMap<Integer, Status> my_status_dict = new HashMap<Integer, Status>();
    HashMap<Integer, Project> my_project_dict = new HashMap<Integer, Project>();
    HashMap<Integer, Chemist> my_chemist_dict = new HashMap<>();
    HashMap<Integer,Vector<Chemist>> chemistsByProject = new HashMap<>();

    Vector<String> super_users = new Vector<String>();

    public static FrontierDAO getInstance(){
        if(_this==null){
            _this = new FrontierDAO();
        }
        return _this;
    }

    private FrontierDAO() {
        try {
            Class.forName("org.postgresql.Driver");
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }
    }

    private void retrieve_My_chemists(){
        my_chemist_dict.clear();
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            statement = connection.prepareStatement("select chemist_id,chemist,project_id from chemists order by chemist_id asc");
            rs = statement.executeQuery();
            while(rs.next()){
                int chemist_id = rs.getInt("chemist_id");
                String chemist = rs.getString("chemist");
                int project_id = rs.getInt("project_id");
                Chemist My_chemist = new Chemist(chemist_id,project_id,chemist);
                my_chemist_dict.put(chemist_id,My_chemist);
                if(!chemistsByProject.containsKey(project_id)){
                    Vector<Chemist> chemists = new Vector<>();
                    chemists.add(new Chemist(0,project_id,""));
                    chemistsByProject.put(project_id, chemists);
                }
                chemistsByProject.get(project_id).add(My_chemist);
            }

        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,rs,statement);
        }
    }

    private void retrieve_status(){
        my_status.clear();
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            statement = connection.prepareStatement("select status_id,status_name from status order by status_id asc");
            rs = statement.executeQuery();
            while(rs.next()){
                int status_id = rs.getInt("status_id");
                String status_name = rs.getString("status_name");
                Status MyStatus = new Status(status_id, status_name);
                my_status.add(MyStatus);
                my_status_dict.put(status_id,MyStatus);
            }

        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,rs,statement);
        }
    }

    private void retrieve_My_projects(){
        my_projects.clear();
        Connection connection = null;
        PreparedStatement statement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            statement = connection.prepareStatement("select project_id,project_name from project order by project_id asc");
            rs = statement.executeQuery();
            while(rs.next()){
                int project_id = rs.getInt("project_id");
                String project_name = rs.getString("project_name");
                Project MyProject = new Project(project_id, project_name);
                my_projects.add(MyProject);
                my_project_dict.put(project_id,MyProject);
            }

        } catch (SQLException e) {
            e.printStackTrace();
        }finally {
            cleanup(connection,rs,statement);
        }

    }

    public Status getStatus(int status_id){
        if(my_status.isEmpty()){
            retrieve_status();
        }
        return my_status_dict.get(status_id);
    }

    public Chemist getAssignedChemist(int chemist_id){
        if(my_chemist_dict.isEmpty()){
            retrieve_My_chemists();
        }
        return my_chemist_dict.get(chemist_id);
    }

    public Vector<Chemist> getChemistsByProject(int project_id){
        if(my_chemist_dict.isEmpty()||chemistsByProject.isEmpty()){
            retrieve_My_chemists();
        }
        return chemistsByProject.get(project_id);
    }

    public Project getMyProject(int project_id){
        if(my_projects.isEmpty()){
            retrieve_My_projects();
        }
        return my_project_dict.get(project_id);
    }

    public Vector<Status> getMy_status() {
        if(my_status.isEmpty()){
            retrieve_status();
        }
        return my_status;
    }

    public Vector<Project> getMy_projects() {
        if(my_projects.isEmpty()){
            retrieve_My_projects();
        }
        return my_projects;
    }

    public Vector<String> getMyProjectNames(){
        if(my_projects.isEmpty()){
            retrieve_My_projects();
        }
        Vector<String> project_names = new Vector<>();
        for(Project pj: my_projects){
            project_names.add(pj.getProject_name());
        }
        return project_names;
    }

/*
    public Vector<String> getSuper_users() {
        if(super_users.isEmpty()){
            Connection connection = null;
            PreparedStatement statement = null;
            ResultSet rs = null;
            try {
                connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary,user,password);
                statement = connection.prepareStatement("select username from super_user");
                rs = statement.executeQuery();
                while(rs.next()){
                    String username = rs.getString("username");
                    super_users.add(username);
                }

            } catch (SQLException e) {
                e.printStackTrace();
            }finally {
                cleanup(connection,rs,statement);
            }
        }
        return super_users;
    }
*/
    public String getDrugName(String smiles) throws SQLException{
        return null;
    }

    public String getCASNumberFromStructure(String smiles) throws SQLException{
        if(smiles==null){
            return null;
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_GoStar, user, password);
            String query = "select a.cas_no from My.cas a inner join rdk.mols b " +
                    "on a.gvk_id = b.gvk_id where b.m<@mol_from_smiles(?::cstring) " +
                    "and b.m@>mol_from_smiles(?::cstring) and b.m::text=?::mol::text;";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smiles);
            preparedStatement.setString(2,smiles);
            preparedStatement.setString(3,smiles);
            rs = preparedStatement.executeQuery();
            String cas_no = null;
            if (rs.next()) {
                cas_no = rs.getString("cas_no");
            }
            if(cas_no!=null){
                return cas_no;
            }
        } finally {
            if(rs!=null){
                rs.close();
            }
            if(preparedStatement!=null){
                preparedStatement.close();
            }
            if(connection!=null){
                connection.close();
            }
        }
        return null;
    }

    public String getSmilesFromDrugName(String drugName) throws SQLException{
        if(drugName.equalsIgnoreCase("name")){
            return null;
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_GoStar, user, password);
            String drugNameUppercase = drugName.toUpperCase();
            String query = "SELECT structure_details.sub_smiles FROM My.compound_synonyms, " +
                           "My.structure_details WHERE structure_details.str_id = compound_synonyms.str_id AND " +
                           "compound_synonyms.synonyms = ? limit 1;";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,drugNameUppercase);
            rs = preparedStatement.executeQuery();
            String smiles = null;
            if (rs.next()) {
                smiles = rs.getString(1);
            }
            if(smiles!=null){
                return smiles;
            }
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
        return null;
    }

    public int getSimilarMoleculeCountFromGVK(String smiles, float similarity, int maxHits) throws SQLException{
        if(maxHits>0){
            return maxHits;
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        int count = 0;
        try {
            connection = DriverManager.getConnection(dbURL_GoStar, user, password);
            String query = "select count(*) from My.structure_details " +
                    "inner join rdk.fps on structure_details.gvk_id = fps.gvk_id " +
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


    public Vector<PropertyMolecule> getSimilarMoleculesFromGVK(String smiles, float similarity, int maxHits) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
        try {
            connection = DriverManager.getConnection(dbURL_GoStar, user, password);
            String query = "select structure_details.sub_smiles,structure_details.gvk_id from My.structure_details " +
                    "inner join rdk.fps on structure_details.gvk_id = fps.gvk_id " +
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
                String sub_smiles = rs.getString("sub_smiles");
                int gvk_id = rs.getInt("gvk_id");
                OEGraphMol mol = ChemFunc.getMolFromMolString(sub_smiles,OEFormat.SMI);
                mol.SetTitle(""+gvk_id);
                oechem.OESetSDData(mol,"gvk_id",""+gvk_id);
                propertyMolecules.add(new PropertyMolecule(mol));
            }
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
        return propertyMolecules;
    }

    public int getSimilarMoleculeCountFromChembl(String smiles, float similarity, int maxHits) throws SQLException{
        if(maxHits>0){
            return maxHits;
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        int count = 0;
        try {
            connection = DriverManager.getConnection(dbURL_chembl, user, password);
            String query = "select count(*) from compound_structures " +
                    "inner join rdk.fps on compound_structures.molregno = fps.molregno " +
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


    public Vector<PropertyMolecule> getSimilarMoleculesFromChembl(String smiles, float similarity, int maxHits) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
        try {
            connection = DriverManager.getConnection(dbURL_chembl, user, password);
            String query = "select compound_structures.canonical_smiles,compound_structures.molregno from compound_structures " +
                    "inner join rdk.fps on compound_structures.molregno = fps.molregno " +
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
                String sub_smiles = rs.getString("canonical_smiles");
                int molregno = rs.getInt("molregno");
                OEGraphMol mol = ChemFunc.getMolFromMolString(sub_smiles,OEFormat.SMI);
                mol.SetTitle(""+molregno);
                oechem.OESetSDData(mol,"molregno",""+molregno);
                propertyMolecules.add(new PropertyMolecule(mol));
            }
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
        return propertyMolecules;
    }

    public int getSimilarMoleculeCountFromAldrich(String smiles, float similarity, int maxHits) throws SQLException{
        if(maxHits>0){
            return maxHits;
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        int count = 0;
        try {
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query = "select count(*) from aldrich_chemicals " +
                    "inner join rdk.fps on aldrich_chemicals.structure_id = fps.structure_id " +
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


    public Vector<PropertyMolecule> getSimilarMoleculesFromAldrich(String smiles, float similarity, int maxHits) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        Vector<PropertyMolecule> propertyMolecules = new Vector<PropertyMolecule>();
        try {
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query = "select aldrich_chemicals.molfile,aldrich_chemicals.structure_id from aldrich_chemicals " +
                    "inner join rdk.fps on aldrich_chemicals.structure_id = fps.structure_id " +
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
                int structure_id = rs.getInt("structure_id");
                OEGraphMol mol = ChemFunc.getMolFromMolString(molfile,OEFormat.SDF);
                mol.SetTitle(""+structure_id);
                oechem.OESetSDData(mol,"structure_id",""+structure_id);
                propertyMolecules.add(new PropertyMolecule(mol));
            }
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
        return propertyMolecules;
    }

    public int getCountFromSmartsFromAldrich(String smarts, int limit) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(molfile) from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and buildingblock = '1' and m@>?::qmol limit %d", limit);
            }else{
                query = "select count(molfile) from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and buildingblock = '1' and m@>?::qmol";
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
            cleanup(connection,rs,preparedStatement);
        }
    }


    public Vector<String> getMolFromSmartsFromAldrich(String smarts, int limit) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query;
            if(limit>0) {
                query = String.format("select molfile,a.structure_id from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and buildingblock = '1' and m@>?::qmol limit %d", limit);
            }else{
                query = "select molfile,a.structure_id from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and buildingblock = '1' and m@>?::qmol";
            }
            System.out.println(query);
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                molFiles.add(molfile);
            }
            return molFiles;
        }catch(Exception e){
            e.printStackTrace();
            return null;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public int getCountFromSmartsFromAldrichScreening(String smarts, int limit) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(molfile) from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and screeningcompound = '1' and m@>?::qmol limit %d", limit);
            }else{
                query = "select count(molfile) from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and screeningcompound = '1' and m@>?::qmol";
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
            cleanup(connection,rs,preparedStatement);
        }
    }


    public Vector<String> getMolFromSmartsFromAldrichScreening(String smarts, int limit) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query;
            if(limit>0) {
                query = String.format("select molfile,a.structure_id from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and screeningcompound = '1' and m@>?::qmol limit %d", limit);
            }else{
                query = "select molfile,a.structure_id from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id and screeningcompound = '1' and m@>?::qmol";
            }
            System.out.println(query);
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                molFiles.add(molfile);
            }
            return molFiles;
        }catch(Exception e){
            e.printStackTrace();
            return null;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }


    public int getCountsFromSmartsFrontier(String smarts,int limit, float clogp, float mw, float psa) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_Frontier, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(molfile) from frontier_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol limit %d", limit);
            }else{
                query = "select count(molfile) from frontier_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol";
            }
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            int count = 0;
            if (rs.next()) {
                count = rs.getInt(1);
            }
            return count;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public Vector<String> getMolFromSmartsFromFrontier(String smarts, int limit, float clogp, float mw, float psa) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL_Frontier, user, password);
            String query;
            if(limit>0) {
                query = String.format("select molfile,acdno,psa,mw,clogp from frontier_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol limit %d", limit);
            }else{
                query = "select molfile,acdno, psa, mw, clogp from frontier_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol";
            }
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                float psa1 = rs.getFloat("psa");
                float mw1 = rs.getFloat("mw");
                float clogp1 = rs.getFloat("clogp");
                OEGraphMol molFromMolString = ChemFunc.getMolFromMolString(molfile, OEFormat.SDF);
                oechem.OESetSDData(molFromMolString,"psa",""+psa1);
                oechem.OESetSDData(molFromMolString,"mw",""+mw1);
                oechem.OESetSDData(molFromMolString,"clogp",""+clogp1);
                molFiles.add(ChemFunc.getMolString(molFromMolString));
            }
            return molFiles;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public int getCountsFromSmartsMarketSelect(String smarts,int limit, float clogp, float mw, float psa) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(molfile) from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol limit %d", limit);
            }else{
                query = "select count(molfile) from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol";
            }
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            int count = 0;
            if (rs.next()) {
                count = rs.getInt(1);
            }
            return count;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public Vector<String> getMolFromSmartsFromMarketSelect(String smarts, int limit, float clogp, float mw, float psa) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL_Aldrich, user, password);
            String query;
            if(limit>0) {
                query = String.format("select molfile,structure_id,psa,mw,clogp from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol limit %d", limit);
            }else{
                query = "select molfile,a.structure_id, psa, mw, clogp from aldrich_chemicals a inner join rdk.mols b on a.structure_id = b.structure_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol";
            }
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                float psa1 = rs.getFloat("psa");
                float mw1 = rs.getFloat("mw");
                float clogp1 = rs.getFloat("clogp");
                OEGraphMol molFromMolString = ChemFunc.getMolFromMolString(molfile, OEFormat.SDF);
                oechem.OESetSDData(molFromMolString,"psa",""+psa1);
                oechem.OESetSDData(molFromMolString,"mw",""+mw1);
                oechem.OESetSDData(molFromMolString,"clogp",""+clogp1);
                molFiles.add(ChemFunc.getMolString(molFromMolString));
            }
            return molFiles;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }


    private String buildPropertyClause(float psa, float mw, float clogp){
        StringBuilder sb = new StringBuilder();
        if(psa>0){
            sb.append(String.format("and psa < %f ",psa));
        }
        if(mw>0){
            sb.append(String.format("and mw < %f ",mw));
        }
        if(clogp>0){
            sb.append(String.format("and clogp < %f ",clogp));
        }
        return sb.toString();
    }

    public int getCountsFromSmartsEnamine(String smarts,int limit, float clogp, float mw, float psa) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_Enamine, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(molfile) from enamine_chemicals a " +
                        "inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        "and m@>?::qmol " +
                        "limit %d", limit);
            }else{
                query = "select count(molfile) from enamine_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol";
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
            cleanup(connection,rs,preparedStatement);
        }
    }

    public Vector<String> getMolFromSmartsFromEnamine(String smarts, int limit, float clogp, float mw, float psa) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL_Enamine, user, password);
            String query;
            if(limit>0) {
                query = String.format("select molfile,catalogue_id,psa,mw,clogp from enamine_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        "and m@>?::qmol " +
                        buildPropertyClause(psa,mw,clogp) +
                        " limit %d", limit);
            }else{
                query = "select molfile,catalogue_id,psa,mw,clogp from enamine_chemicals a inner join rdk.mols b " +
                        "on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        "and m@>?::qmol";
            }
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                float psa1 = rs.getFloat("psa");
                float mw1 = rs.getFloat("mw");
                float clogp1 = rs.getFloat("clogp");
                OEGraphMol molFromMolString = ChemFunc.getMolFromMolString(molfile, OEFormat.SDF);
                oechem.OESetSDData(molFromMolString,"psa",""+psa1);
                oechem.OESetSDData(molFromMolString,"mw",""+mw1);
                oechem.OESetSDData(molFromMolString,"clogp",""+clogp1);
                molFiles.add(ChemFunc.getMolString(molFromMolString));
            }
            return molFiles;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public int getCountsFromSmartsPharmaron(String smarts,int limit, float clogp, float mw, float psa) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(dbURL_Pharmaron, user, password);
            String query;
            if(limit>0) {
                query = String.format("select count(molfile) from enamine_chemicals a " +
                        "inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        "and m@>?::qmol " +
                        "limit %d", limit);
            }else{
                query = "select count(molfile) from pharmaron_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        " and m@>?::qmol";
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
            cleanup(connection,rs,preparedStatement);
        }
    }

    public Vector<String> getMolFromSmartsFromPharmaron(String smarts, int limit, float clogp, float mw, float psa) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            Vector<String> molFiles = new Vector<String>();
            connection = DriverManager.getConnection(dbURL_Pharmaron, user, password);
            String query;
            if(limit>0) {
                query = String.format("select molfile,catalogue_id,psa,mw,clogp from pharmaron_chemicals a inner join rdk.mols b on a.chemical_id = b.chemical_id " +
                        "and m@>?::qmol " +
                        buildPropertyClause(psa,mw,clogp) +
                        " limit %d", limit);
            }else{
                query = "select molfile,catalogue_id,psa,mw,clogp from pharmaron_chemicals a inner join rdk.mols b " +
                        "on a.chemical_id = b.chemical_id " +
                        buildPropertyClause(psa,mw,clogp) +
                        "and m@>?::qmol";
            }
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while (rs.next()) {
                String molfile = rs.getString("molfile");
                float psa1 = rs.getFloat("psa");
                float mw1 = rs.getFloat("mw");
                float clogp1 = rs.getFloat("clogp");
                OEGraphMol molFromMolString = ChemFunc.getMolFromMolString(molfile, OEFormat.SDF);
                oechem.OESetSDData(molFromMolString,"psa",""+psa1);
                oechem.OESetSDData(molFromMolString,"mw",""+mw1);
                oechem.OESetSDData(molFromMolString,"clogp",""+clogp1);
                molFiles.add(ChemFunc.getMolString(molFromMolString));
            }
            return molFiles;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }


    public Vector<CoreStatus> getAllCoreStatus() throws SQLException{
        Vector<CoreStatus> all_status = new Vector<CoreStatus>();
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try{
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            preparedStatement = connection.prepareStatement("select status_id,status_name from status order by status_id asc");
            rs = preparedStatement.executeQuery();
            while(rs.next()){
                if(!rs.wasNull()){
                    CoreStatus status = new CoreStatus(rs.getInt("status_id"), rs.getString("status_name"));
                    all_status.add(status);
                }
            }
            return  all_status;
        }finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public String getCoreStatus(int status_id) throws SQLException{
        if(!statusDict.containsKey(status_id)) {
            Connection connection = null;
            PreparedStatement preparedStatement = null;
            ResultSet rs = null;
            try {
                connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
                preparedStatement = connection.prepareStatement("select status_id, status_name from status");
                rs = preparedStatement.executeQuery();
                while (rs.next()) {
                    if (!rs.wasNull()) {
                        int status_id1 = rs.getInt("status_id");
                        String status_name = rs.getString("status_name");
                        statusDict.put(status_id1,status_name);
                    }
                }
                if(statusDict.containsKey(status_id)) {
                    return statusDict.get(status_id);
                }else{
                    return "";
                }
            } finally {
                cleanup(connection, rs, preparedStatement);
            }
        }else{
            return statusDict.get(status_id);
        }
    }

    public void insertSampleMolecules(int core_id, Vector<PropertyMolecule> propertyMolecules, ProgressReporter progressReporter) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            if(propertyMolecules==null||propertyMolecules.isEmpty()){
                throw new SQLException("No molecules available.");
            }
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            connection.setAutoCommit(false);
            String query = "insert into core_enumeration (core_id, compound_sdf, unique_smiles) values (?,?,?)";
            preparedStatement = connection.prepareStatement(query);
            int n=0;
            for(PropertyMolecule pmol:propertyMolecules) {
                if(progressReporter!=null){
                    progressReporter.reportProgress("Inserting molecules ...",100*n/propertyMolecules.size());
                }
                String molStr = OEChemFunc.getInstance().getStringFromOEMol(pmol.getMol(), OEFormat.SDF);
                String smiles = oechem.OEMolToSmiles(pmol.getMol());
                preparedStatement.setInt(1, core_id);
                preparedStatement.setString(2,molStr);
                preparedStatement.setString(3,smiles);
                preparedStatement.addBatch();
                n+=1;
            }
            preparedStatement.executeBatch();
            connection.commit();
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public Vector<PropertyMolecule> getSampleMolecules(int core_id) throws SQLException{
        Vector<PropertyMolecule> molecules = new Vector<PropertyMolecule>();
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try{
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query = "select compound_sdf from core_enumeration where core_id = ?";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,core_id);
            rs = preparedStatement.executeQuery();
            while(rs.next()){
                String compoundSdf = rs.getString("compound_sdf");
                OEGraphMol mol = OEChemFunc.getInstance().convertString2OEMol(compoundSdf,OEFormat.MDL);
                molecules.add(new PropertyMolecule(mol));
            }
            ChemFunc.calculateOEProperty(molecules);
        }finally{
            cleanup(connection,rs,preparedStatement);
        }
        return molecules;
    }

    public void removeMyCore(int core_id) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query = "delete from cores where core_id = ?";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,core_id);
            preparedStatement.executeUpdate();
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public void insertOrUpdateMyComoundRanking(int ranking, int mole_id, String username) throws SQLException{
        int old_ranking = getMyCompoundRanking(mole_id, username);
        if(old_ranking==-1){
            insertMyComoundRanking(ranking,mole_id,username);
        }else if(old_ranking!=ranking){
            updateMyCompoundRanking(ranking,mole_id,username);
        }
    }

    public void deleteCompoundRanking(int mole_id, String username) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_CompoundTracker, user, password);
            String query = "delete from rankings where mole_id = ? and username = ?";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,mole_id);
            preparedStatement.setString(2,username);
            preparedStatement.executeUpdate();
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }



    private int getMyCompoundRanking(int mole_id, String username) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            preparedStatement = connection.prepareStatement("select ranking from rankings where mole_id = ? and username = ? ");
            preparedStatement.setInt(1,mole_id);
            preparedStatement.setString(2,username);
            rs = preparedStatement.executeQuery();
            if(rs.next()) {
                int ranking = rs.getInt("ranking");
                if (rs.wasNull()) {
                    return -1;
                } else {
                    return ranking;
                }
            }else{
                return -1;
            }
        }finally{
            cleanup(connection, rs,preparedStatement);
        }
    }

    private void insertMyComoundRanking(int ranking, int mole_id, String username ) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            preparedStatement = connection.prepareStatement("insert into rankings (mole_id,ranking,username) " +
                    "values (?,?,?)");
            preparedStatement.setInt(1,mole_id);
            preparedStatement.setInt(2,ranking);
            preparedStatement.setString(3,username);
            preparedStatement.executeUpdate();
        }finally{
            cleanup(connection,null,preparedStatement);
        }
    }

    private void updateMyCompoundRanking(int ranking, int mole_id, String username) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            preparedStatement = connection.prepareStatement("update rankings set ranking = ? where mole_id = ? and username = ? ");
            preparedStatement.setInt(1,ranking);
            preparedStatement.setInt(2,mole_id);
            preparedStatement.setString(3,username);
            preparedStatement.executeUpdate();
        }finally{
            cleanup(connection,null,preparedStatement);
        }
    }

    public void deleteMyCompounds(Vector<Compound> MyCompounds) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            preparedStatement = connection.prepareStatement(
                    "delete from molecules where mole_id = ?");
            connection.setAutoCommit(false);
            for(Compound c:MyCompounds){
                preparedStatement.setInt(1,c.getId());
                preparedStatement.addBatch();
            }
            preparedStatement.executeBatch();
            connection.commit();
        }finally{
            cleanup(connection,null,preparedStatement);
        }
    }

    public void updateMyCompound(Compound MyCompound) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            preparedStatement = connection.prepareStatement("update molecules set molfile=?, mol_name=?, status_id=?, document=?, smiles=?, chemist_id=?, psa=?, mw=?, clogp=?, cns_mpo=?, cns_tempo=? where mole_id = ?");
            int mol_id = MyCompound.getId();
            String mol_name = MyCompound.getName();
            OEGraphMol mol = MyCompound.getPropertyMol().getMol();
            Molecule chemaxonMol = OEChemFunc.getInstance().convertOEChemMol(mol);
            String molStr = ChemFunc.getMolString(mol);
            String smiles;
            try {
                smiles = MolExporter.exportToFormat(chemaxonMol,"smiles:u,a");
            } catch (IOException e) {
                e.printStackTrace();
                throw new SQLException(e.getMessage());
            }
            byte[] cdx = new byte[0];
            try {
                cdx = MolExporter.exportToBinFormat(chemaxonMol, "cdx");
            } catch (IOException e) {
                e.printStackTrace();
                throw new SQLException(e.getMessage());
            }
            preparedStatement.setString(1,molStr);
            preparedStatement.setString(2,mol_name);
            preparedStatement.setInt(3,MyCompound.getStatus_id());
            preparedStatement.setBinaryStream(4,new ByteArrayInputStream(cdx),cdx.length);
            preparedStatement.setString(5, smiles);
            Chemist assignedChemist = MyCompound.getAssignedChemist();
            if(assignedChemist==null){
                preparedStatement.setNull(6, Types.INTEGER);
            }else {
                preparedStatement.setInt(6, assignedChemist.getChemist_id());
            }
            preparedStatement.setDouble(7,MyCompound.getPropertyMol().getPSA());
            preparedStatement.setDouble(8,MyCompound.getPropertyMol().getMW());
            preparedStatement.setDouble(9,MyCompound.getPropertyMol().getCLogP());
            MolProperty cns_mpo = MyCompound.getPropertyMol().getProperty("CNS MPO");
            preparedStatement.setDouble(10, cns_mpo==null?0.0:cns_mpo.getValue());
            MolProperty cns_mTEMPO = MyCompound.getPropertyMol().getProperty("CNS mTEMPO");
            preparedStatement.setDouble(11, cns_mTEMPO==null?0:cns_mTEMPO.getValue());

            preparedStatement.setInt(12,mol_id);
            preparedStatement.executeUpdate();
        }finally{
            cleanup(connection,null,preparedStatement);
        }
        if(MyCompound.isRankingChanged()) {
            insertOrUpdateMyComoundRanking(MyCompound.getRank(), MyCompound.getId(), InSlilicoPanel.getInstance().getUserName());
        }
    }

    public int insertMyCompound(Compound MyCompound) throws SQLException{
        int regid = -1;
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            preparedStatement = connection.prepareStatement(
                    "insert into molecules (molfile, mol_name, project_id,status_id,chemist,document, smiles, chemist_id, psa, mw, clogp, cns_mpo, cns_tempo) " +
                            "(select ?,?,?,?,?,?,?,?,?,?,?,?,? where not exists (select smiles,project_id from molecules where smiles = ? and project_id = ?) )", PreparedStatement.RETURN_GENERATED_KEYS);
            String mol_name = MyCompound.getName();
            OEGraphMol mol = MyCompound.getPropertyMol().getMol();
            Molecule chemaxonMol = OEChemFunc.getInstance().convertOEChemMol(mol);
            String molStr = ChemFunc.getMolString(mol);
            String smiles;
            try {
                smiles = MolExporter.exportToFormat(chemaxonMol,"smiles:u,a");
            } catch (IOException e) {
                e.printStackTrace();
                return -1;
            }
            byte[] cdx = new byte[0];
            try {
                cdx = MolExporter.exportToBinFormat(chemaxonMol, "cdx");
            } catch (IOException e) {
                e.printStackTrace();
                return -1;
            }
            preparedStatement.setString(1,molStr);
            preparedStatement.setString(2,mol_name);
            preparedStatement.setInt(3,MyCompound.getProject_id());
            preparedStatement.setInt(4,MyCompound.getStatus_id());
            preparedStatement.setString(5,MyCompound.getChemist());
            preparedStatement.setBinaryStream(6,new ByteArrayInputStream(cdx),cdx.length);
            preparedStatement.setString(7, smiles);
            Chemist assignedChemist = MyCompound.getAssignedChemist();
            if(assignedChemist==null){
                preparedStatement.setNull(8, Types.INTEGER);
            }else {
                preparedStatement.setInt(8, assignedChemist.getChemist_id());
            }
            preparedStatement.setDouble(9,MyCompound.getPropertyMol().getPSA());
            preparedStatement.setDouble(10,MyCompound.getPropertyMol().getMW());
            preparedStatement.setDouble(11,MyCompound.getPropertyMol().getCLogP());
            preparedStatement.setDouble(12,MyCompound.getPropertyMol().getProperty("CNS MPO").getValue());
            preparedStatement.setDouble(13,MyCompound.getPropertyMol().getProperty("CNS mTEMPO").getValue());
            preparedStatement.setString(14,smiles);
            preparedStatement.setInt(15,MyCompound.getProject_id());
            preparedStatement.executeUpdate();
            rs = preparedStatement.getGeneratedKeys();
            if (rs.next()) {
                int id = rs.getInt(1);
                if(!rs.wasNull()) {
                    regid = id;
                }
            }
            return regid;

        }finally{
            cleanup(connection,rs,preparedStatement);
        }
    }

    public Vector<Integer> insertMyCompoundsBatch(String sdfName, String chemist, String nameTag, int project_id, int status_id, Chemist MyChemist) throws SQLException{
        if(Strings.isNullOrEmpty(sdfName)){
            throw new SQLException("No SDF defined.");
        }
        if(Strings.isNullOrEmpty(chemist)){
            throw new SQLException("No chemist defined.");
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try{
            Vector<Integer> idList = new Vector<Integer>();
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            connection.setAutoCommit(false);
            preparedStatement = connection.prepareStatement(
                    "insert into molecules (molfile, mol_name, project_id,status_id,chemist,document, smiles, chemist_id) " +
                            "(select ?,?,?,?,?,?,?,? where not exists (select smiles,project_id from molecules where smiles = ? and project_id = ?) )", PreparedStatement.RETURN_GENERATED_KEYS);
            oemolistream ifs = new oemolistream();
            ifs.open(sdfName);
            OEGraphMol mol = new OEGraphMol();
            while(oechem.OEReadMolecule(ifs,mol)){
                String mol_name = "";
                if(nameTag!=null&&(!nameTag.trim().isEmpty())&&oechem.OEHasSDData(mol,nameTag)){
                    mol_name = oechem.OEGetSDData(mol,nameTag);
                }else{
                    mol_name = mol.GetTitle();
                    if(mol_name.isEmpty()) {
                        mol_name = "Anonymous";
                    }
                }
                Molecule chemaxonMol = OEChemFunc.getInstance().convertOEChemMol(mol);
                String molStr = ChemFunc.getMolString(mol);
                String smiles;
                try {
                    smiles = MolExporter.exportToFormat(chemaxonMol,"smiles:u,a");
                } catch (IOException e) {
                    e.printStackTrace();
                    continue;
                }
                byte[] cdx = new byte[0];
                try {
                    cdx = MolExporter.exportToBinFormat(chemaxonMol, "cdx");
                } catch (IOException e) {
                    e.printStackTrace();
                }
                if(cdx==null){
                    continue;
                }
//molfile,project_id,status_id,chemist,document, smiles
                preparedStatement.setString(1,molStr);
                preparedStatement.setString(2,mol_name);
                preparedStatement.setInt(3,project_id);
                preparedStatement.setInt(4,status_id);
                preparedStatement.setString(5,chemist);
                preparedStatement.setBinaryStream(6,new ByteArrayInputStream(cdx),cdx.length);
                preparedStatement.setString(7, smiles);
                if(MyChemist==null||MyChemist.getChemist_id()==0){
                    preparedStatement.setNull(8, Types.INTEGER);
                }else {
                    preparedStatement.setInt(8, MyChemist.getChemist_id());
                }

                preparedStatement.setString(9,smiles);
                preparedStatement.setInt(10,project_id);
                preparedStatement.addBatch();
            }
            preparedStatement.executeBatch();
            connection.commit();
            rs = preparedStatement.getGeneratedKeys();
            while (rs.next()) {
                int id = rs.getInt(1);
                if(!rs.wasNull()) {
                    idList.add(id);
                }
            }
            return idList;

        }finally{
            cleanup(connection,rs,preparedStatement);
        }

    }

    public void addComment(int mol_id, String chemist, String comment) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            String query = "insert into comments (molecule_id, chemist, comment, date) values (?,?,?,now())";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,mol_id);
            preparedStatement.setString(2,chemist);
            preparedStatement.setString(3,comment);
            preparedStatement.executeUpdate();
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }

    public void deleteComment(int mol_id, String chemist) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try{
            connection = DriverManager.getConnection(getDbURL_CompoundTracker,user,password);
            String query = "delete from comments where molecule_id = ? and chemist = ?";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,mol_id);
            preparedStatement.setString(2,chemist);
            preparedStatement.executeUpdate();
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }

    public String getMyCompoundComment(int mol_id) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_CompoundTracker, user, password);
            String query = "select chemist,comment from comments where molecule_id = ? order by date asc";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,mol_id);
            rs = preparedStatement.executeQuery();
            StringBuilder sb = new StringBuilder();
            while(rs.next()){
                String chemist = rs.getString("chemist");
                String comment = rs.getString("comment");
                sb.append(String.format("%s:%s\n",chemist.trim(),comment.trim()));
            }
            return sb.toString();
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }


    public Vector<Compound> getMyCompounds(int project_id) throws SQLException{
        Vector<Compound> MyCompounds = new Vector<Compound>();
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_CompoundTracker, user, password);
            String query = "select m.mole_id,m.mol_name, m.document, m.molfile,m.status_id,m.smiles, m.document," +
                    "m.date, m.chemist, m.chemist_id, avg(r.ranking) as rankings, m.psa, m.mw, m.clogp, m.cns_mpo, m.cns_tempo from molecules m " +
                    "full join rankings r on m.mole_id = r.mole_id "+
                    "where m.project_id = ? group by m.mole_id "+
                    "order by date desc";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,project_id);
            rs = preparedStatement.executeQuery();
            while(rs.next()){
                int mole_id = rs.getInt("mole_id");
                String name = rs.getString("mol_name");
                String smiles = rs.getString("smiles");
                String molecule = rs.getString("molfile");
                byte[] document = rs.getBytes("document");
                String chemist = rs.getString("chemist");
                Date date = rs.getDate("date");
                int status_id = rs.getInt("status_id");
                int rank = rs.getInt("rankings");
                int chemist_id  = rs.getInt("chemist_id");
                double psa = rs.getDouble("psa");
                double mw = rs.getDouble("mw");
                double clogp = rs.getDouble("clogp");
                double cns_mpo = rs.getDouble("cns_mpo");
                double cns_tempo = rs.getDouble("cns_tempo");
                Compound compound = new Compound(mole_id,name,molecule,smiles,document,chemist,date,project_id,status_id,rank,chemist_id,clogp,mw,psa,cns_tempo,cns_mpo);
                MyCompounds.add(compound);
            }
            return MyCompounds;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }


    public Vector<Core> getMyCoresBySubstructure(String smarts) throws SQLException{
        Vector<Core> cores = new Vector<Core>();
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query = "select * from cores inner join status on cores.status_id = status.status_id where core_m@>?::qmol and active = 1 order by date desc";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,smarts);
            rs = preparedStatement.executeQuery();
            while(rs.next()){
                int core_id = rs.getInt("core_id");
                String coreMol = rs.getString("core_mol");
                String name = rs.getString("name");
                String coreSmiles = rs.getString("core_smiles");
                int num_of_diversity = rs.getInt("num_of_diversity");
                String description = rs.getString("description");
                byte[] document = rs.getBytes("document");
                String chemist = rs.getString("chemist");
                Date date = rs.getDate("date");
                String source = rs.getString("source");
                String comment = rs.getString("comment");
                int status_id = rs.getInt("status_id");
                String status_name = rs.getString("status_name");
                Core core = new Core(core_id,name, coreMol,coreSmiles,description,document,chemist,date,source,comment,status_id);
                cores.add(core);
            }
            return cores;
        } finally {
            cleanup(connection,rs,preparedStatement);
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



    public Vector<Core> getAllMyCores() throws SQLException{
        Vector<Core> cores = new Vector<Core>();
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query = "select * from cores inner join status on cores.status_id = status.status_id where active = 1 order by date desc";
            preparedStatement = connection.prepareStatement(query);
            rs = preparedStatement.executeQuery();
            while(rs.next()){
                int core_id = rs.getInt("core_id");
                String name = rs.getString("name");
                String coreMol = rs.getString("core_mol");
                String coreSmiles = rs.getString("core_smiles");
                int num_of_diversity = rs.getInt("num_of_diversity");
                String description = rs.getString("description");
                byte[] document = rs.getBytes("document");
                String chemist = rs.getString("chemist");
                Date date = rs.getDate("date");
                String source = rs.getString("source");
                String comment = rs.getString("comment");
                int status_id = rs.getInt("status_id");
                Core core = new Core(core_id,name, coreMol,coreSmiles,description,document,chemist,date,source,comment,status_id);
                cores.add(core);
            }
            return cores;
        } finally {
            cleanup(connection,rs,preparedStatement);
        }
    }

    public void insertMyCoresBatchWrapper(String sdfPathName, String chemist, String source, String libraryNameTag) throws SQLException {
            Vector<Integer> idList = insertMyCoresBatch(sdfPathName, chemist, source, libraryNameTag);
//            insertCoreToMolBatch(idList);
//            insertFPforCoreBatch(idList);
    }


    private Vector<Integer> insertMyCoresBatch(String sdfName, String chemist, String source, String tagName) throws SQLException{
        if(Strings.isNullOrEmpty(sdfName)){
            throw new SQLException("No SDF defined.");
        }
        if(Strings.isNullOrEmpty(chemist)){
            throw new SQLException("No chemist defined.");
        }
        if(Strings.isNullOrEmpty(source)){
            source = "Cellarity";
        }
        if(Strings.isNullOrEmpty(tagName)){
            throw new SQLException("No name tag defined. ");
        }
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs = null;
        try{
            Vector<Integer> idList = new Vector<Integer>();
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary,user,password);
            connection.setAutoCommit(false);
            preparedStatement = connection.prepareStatement("insert into cores (name, core_mol,core_smiles,num_of_diversity,description,document,chemist, " +
                    "source, core_m, torsionbv, mfp2, ffp2 " +
                    "values (?,?,?,?,?,?,?,?,mol_from_smiles(?::cstring),torsionbv_fp(mol_from_smiles(?::cstring))," +
                    "morganbv_fp(mol_from_smiles(?::cstring)),featmorganbv_fp(mol_from_smiles(?::cstring)))", PreparedStatement.RETURN_GENERATED_KEYS);
            oemolistream ifs = new oemolistream();
            ifs.open(sdfName);
            OEGraphMol mol = new OEGraphMol();
            while(oechem.OEReadMolecule(ifs,mol)){
                if(oechem.OEHasSDData(mol,tagName)){
                    int number_of_diversities = 0;
                    for(OEAtomBase atm:mol.GetAtoms()){
                        if(new OEIsRGroup(0).constCall(atm)){
                            number_of_diversities += 1;
                        }
                    }

                    Molecule chemaxonMol = OEChemFunc.getInstance().convertOEChemMol(mol);
                    String molStr = ChemFunc.getMolString(mol);
                    String coreSmiles;
                    try {
                        coreSmiles = MolExporter.exportToFormat(chemaxonMol,"smiles:u,a");
                    } catch (IOException e) {
                        e.printStackTrace();
                        continue;
                    }
                    byte[] cdx = new byte[0];
                    try {
                        cdx = MolExporter.exportToBinFormat(chemaxonMol, "cdx");
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    if(cdx==null){
                        continue;
                    }

                    String name = oechem.OEGetSDData(mol,tagName);
                    preparedStatement.setString(1,name);
                    preparedStatement.setString(2,molStr);
                    preparedStatement.setString(3,coreSmiles);
                    preparedStatement.setInt(4,number_of_diversities);
                    preparedStatement.setString(5,"Automatically generated text.");
                    preparedStatement.setBinaryStream(6,new ByteArrayInputStream(cdx),cdx.length);
                    preparedStatement.setString(7, chemist);
                    preparedStatement.setString(8,source);
                    preparedStatement.setString(9,coreSmiles);
                    preparedStatement.setString(10,coreSmiles);
                    preparedStatement.setString(11,coreSmiles);
                    preparedStatement.setString(12,coreSmiles);
                    preparedStatement.addBatch();
                }else{
                    System.err.println(String.format("Failed to find tag in molecule %s, skipped.",mol.GetTitle()));
                }
            }
            preparedStatement.executeBatch();
            connection.commit();
            rs = preparedStatement.getGeneratedKeys();
            while (rs.next()) {
                int id = rs.getInt(1);
                if(!rs.wasNull()) {
                    idList.add(id);
                }
            }
            return idList;

        }finally{
            cleanup(connection,rs,preparedStatement);
        }

    }
    /*
    private void insertCoreToMolBatch(Vector<Integer> core_ids) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            connection.setAutoCommit(false);
            String query = "insert into rdk.mols (core_id,m)  values (?,mol_from_smiles((select core_smiles from cores where core_id = ?)::cstring))";
            preparedStatement = connection.prepareStatement(query);
            for(int core_id:core_ids) {
                preparedStatement.setInt(1, core_id);
                preparedStatement.setInt(2, core_id);
                preparedStatement.addBatch();
            }
            preparedStatement.executeBatch();
            connection.commit();
        }finally {
            cleanup(connection,null,preparedStatement);
        }

    }


    private void insertFPforCoreBatch(Vector<Integer> core_ids) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            connection.setAutoCommit(false);
            String query;
            query = "insert into rdk.fps (core_id, torsionbv, mfp2, ffp2) values (?," +
                    "torsionbv_fp(mol_from_smiles((select core_smiles from cores where core_id = ?)::cstring)), " +
                    "morganbv_fp(mol_from_smiles((select core_smiles from cores where core_id = ?)::cstring)), " +
                    "featmorganbv_fp(mol_from_smiles((select core_smiles from cores where core_id = ?)::cstring)))";
            preparedStatement = connection.prepareStatement(query);
            for(int core_id:core_ids) {
                preparedStatement.setInt(1, core_id);
                preparedStatement.setInt(2, core_id);
                preparedStatement.setInt(3, core_id);
                preparedStatement.setInt(4, core_id);
                preparedStatement.addBatch();
            }
            preparedStatement.executeBatch();
            connection.commit();
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }
    */

    public void updateMyCoreWrapper(int core_id, String name, String coreMol, String coreSmiles, int numOfDiversity, String description, byte[] cdx, String chemist, String source, int status_id) throws SQLException{
        if(Strings.isNullOrEmpty(coreMol)||Strings.isNullOrEmpty(coreSmiles)){
            throw new SQLException("No core defined.");
        }
        if(numOfDiversity==0){
            throw new SQLException("No R group defined.");
        }
        if(cdx==null||cdx.length==0){
            throw new SQLException("No reaction defined.");
        }
        if(Strings.isNullOrEmpty(chemist)){
            chemist = "Anonymous";
        }
        updateMyCore(core_id, name,coreMol, coreSmiles, numOfDiversity, description, cdx, chemist, source, status_id);
//        updateCoreToMol(core_id,coreSmiles);
//        updateFPforCore(core_id,coreSmiles);
    }


    private void updateMyCore(int core_id, String name, String coreMol, String coreSmiles, int numOfDiversity, String description, byte[] cdx, String chemist, String source, int status_id) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query;
            //torsionbv_fp(mol_from_smiles(core_smiles)), morganbv_fp(mol_from_smiles((core_smiles)), " +
            //                    "featmorganbv_fp(mol_from_smiles((core_smiles))
            query = "update cores set name = ?, core_mol = ?,core_smiles = ?,num_of_diversity = ?,description = ?," +
                    "document = ?,chemist = ?, source = ?, status_id = ?, " +
                    "core_m = mol_from_smiles(?::cstring), torsionbv = torsionbv_fp(mol_from_smiles(?::cstring))," +
                    " mfp2 = morganbv_fp(mol_from_smiles(?::cstring)), ffp2 = featmorganbv_fp(mol_from_smiles(?::cstring)) where core_id = ?";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setString(1,name);
            preparedStatement.setString(2,coreMol);
            preparedStatement.setString(3,coreSmiles);
            preparedStatement.setInt(4,numOfDiversity);
            preparedStatement.setString(5,description);
            preparedStatement.setBinaryStream(6, new ByteArrayInputStream(cdx),cdx.length);
            preparedStatement.setString(7, chemist);
            preparedStatement.setString(8,source);
            preparedStatement.setInt(9,status_id);
            preparedStatement.setString(10,coreSmiles);
            preparedStatement.setString(11,coreSmiles);
            preparedStatement.setString(12,coreSmiles);
            preparedStatement.setString(13,coreSmiles);
            preparedStatement.setInt(14,core_id);
            preparedStatement.executeUpdate();
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }

//    private void updateCoreToMol(int core_id, String coreSmiles) throws SQLException{
//        Connection connection = null;
//        PreparedStatement preparedStatement = null;
//        try {
//            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
//            String query;
//            query = "update rdk.mols set m = mol_from_smiles(?::cstring) where core_id = ?";
//            preparedStatement = connection.prepareStatement(query);
//            preparedStatement.setString(1,coreSmiles);
//            preparedStatement.setInt(2,core_id);
//            preparedStatement.executeUpdate();
//        }finally {
//            cleanup(connection,null,preparedStatement);
//        }
//    }

//    private void updateFPforCore(int core_id, String coreSmiles) throws SQLException{
//        Connection connection = null;
//        PreparedStatement preparedStatement = null;
//        try {
//            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
//            String query;
//            query = "update rdk.fps set torsionbv = torsionbv_fp(mol_from_smiles(?::cstring)), " +
//                    "mfp2 = morganbv_fp(mol_from_smiles(?::cstring)), " +
//                    "ffp2 = featmorganbv_fp(mol_from_smiles(?::cstring)) where core_id = ?";
//            preparedStatement = connection.prepareStatement(query);
//            preparedStatement.setString(1,coreSmiles);
//            preparedStatement.setString(2,coreSmiles);
//            preparedStatement.setString(3,coreSmiles);
//            preparedStatement.setInt(4,core_id);
//            preparedStatement.executeUpdate();
//        }finally {
//            cleanup(connection,null,preparedStatement);
//        }
//    }

    public Core getCore(int core_id) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        ResultSet rs;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            preparedStatement = connection.prepareStatement("select * from cores where core_id = ?");
            preparedStatement.setInt(1,core_id);
            rs = preparedStatement.executeQuery();
            if(rs.next()){
                //name, core_mol,core_smiles,num_of_diversity,description,document,chemist, source
                String name = rs.getString("name");
                String coreMol = rs.getString("core_mol");
                String coreSmiles = rs.getString("core_smiles");
                int num_of_diversity = rs.getInt("num_of_diversity");
                String description = rs.getString("description");
                byte[] document = rs.getBytes("document");
                String chemist = rs.getString("chemist");
                Date date = rs.getDate("date");
                String source = rs.getString("source");
                String comment = rs.getString("comment");
                int status_id = rs.getInt("status_id");
                return new Core(core_id,name, coreMol,coreSmiles,description,document,chemist,date,source,comment,status_id);
            }else{
                return null;
            }
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }


    public int insertMyCoreWrapper(String name, String coreMol, String coreSmiles, int numOfDiversity, String description, byte[] cdx, String chemist, String source,int status_id) throws SQLException{
        if(Strings.isNullOrEmpty(coreMol)||Strings.isNullOrEmpty(coreSmiles)){
            throw new SQLException("No core defined.");
        }
        if(numOfDiversity==0){
            throw new SQLException("No R group defined.");
        }
        if(cdx==null||cdx.length==0){
            throw new SQLException("No reaction defined.");
        }
        if(Strings.isNullOrEmpty(chemist)){
            chemist = "Anonymous";
        }
        int id = insertMyCore(name, coreMol,coreSmiles,numOfDiversity,description,cdx,chemist,source, status_id);
        return id;
    }

    private int insertMyCore(String name, String coreMol, String coreSmiles, int numOfDiversity, String description, byte[] cdx, String chemist, String source, int status_id) throws SQLException {
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query;
            query = "insert into cores (name, core_mol,core_smiles,num_of_diversity,description,document,chemist, source, status_id,core_m,torsionbv,mfp2,ffp2)" +
                    " values (?,?,?,?,?,?,?,?,?,mol_from_smiles(?::cstring), torsionbv_fp(mol_from_smiles(?::cstring)), " +
                    "morganbv_fp(mol_from_smiles(?::cstring)), featmorganbv_fp(mol_from_smiles(?::cstring)))";
            preparedStatement = connection.prepareStatement(query,Statement.RETURN_GENERATED_KEYS);
            preparedStatement.setString(1,name);
            preparedStatement.setString(2,coreMol);
            preparedStatement.setString(3,coreSmiles);
            preparedStatement.setInt(4,numOfDiversity);
            preparedStatement.setString(5,description);
            preparedStatement.setBinaryStream(6, new ByteArrayInputStream(cdx),cdx.length);
            preparedStatement.setString(7, chemist);
            preparedStatement.setString(8,source);
            preparedStatement.setInt(9,status_id);
            preparedStatement.setString(10,coreSmiles);
            preparedStatement.setString(11,coreSmiles);
            preparedStatement.setString(12,coreSmiles);
            preparedStatement.setString(13,coreSmiles);
            int result = preparedStatement.executeUpdate();
            if(result == 1){
                ResultSet generatedKeys = preparedStatement.getGeneratedKeys();
                if (generatedKeys.next()) {
                    return generatedKeys.getInt(1);
                }
                else {
                    throw new SQLException("Creating core failed, no ID obtained.");
                }
            }else{
                throw new SQLException("Creating core failed, no ID obtained.");
            }
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }

//            "insert into rdk.mols (bio_number,m) select %s, mol_from_ctab(%s) where not exists (select bio_number from rdk.mols where bio_number = %s)"
//            "insert into rdk.fps (bio_number, torsionbv, mfp2, ffp2) select %s, torsionbv_fp(mol_from_ctab(%s)), morganbv_fp(mol_from_ctab(%s)), featmorganbv_fp(mol_from_ctab(%s)) where not exists (select bio_number from rdk.fps where bio_number = %s)"

    private void insertCoreToMol(int core_id, String coreSmiles) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query;
            query = "insert into rdk.mols (core_id,m)  values (?,mol_from_smiles(?::cstring))";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,core_id);
            preparedStatement.setString(2,coreSmiles);
            preparedStatement.executeUpdate();
        }finally {
            cleanup(connection,null,preparedStatement);
        }

    }

    private void insertFPforCore(int core_id, String coreSmiles) throws SQLException{
        Connection connection = null;
        PreparedStatement preparedStatement = null;
        try {
            connection = DriverManager.getConnection(getDbURL_VirtualCmpdLibrary, user, password);
            String query;
            query = "insert into rdk.fps (core_id, torsionbv, mfp2, ffp2) values (?,torsionbv_fp(mol_from_smiles(?::cstring)), " +
                    "morganbv_fp(mol_from_smiles(?::cstring)), " +
                    "featmorganbv_fp(mol_from_smiles(?::cstring)))";
            preparedStatement = connection.prepareStatement(query);
            preparedStatement.setInt(1,core_id);
            preparedStatement.setString(2,coreSmiles);
            preparedStatement.setString(3,coreSmiles);
            preparedStatement.setString(4,coreSmiles);
            preparedStatement.executeUpdate();
        }finally {
            cleanup(connection,null,preparedStatement);
        }
    }






    public static void main(String[] args) {
        try {

            FrontierDAO.getInstance().insertOrUpdateMyComoundRanking(4,3235,"jfeng1");
            //FrontierDAO.getInstance().deleteCompoundRanking(3235,"jfeng1");
//            int count = FrontierDAO.getInstance().getSimilarMoleculeCountFromChembl("Clc1nc(cnc1)NC2CCN(CC2)Cc3sc(nc3)C",0.4F,-1);
//            String comment = FrontierDAO.getInstance().getMyCompoundComment(3236);
//            FrontierDAO.getInstance().addComment(3236,"jfeng1","Very Good!");
//            System.out.println(comment);

        } catch (SQLException e) {
            e.printStackTrace();
        }
//        String errorMsg = "dfaetse \n Detail: Key (core_smiles)=([*]c1ccc([*])cc1) already exists.\n";
//        String errorMsg2 = errorMsg.replace("\n","");
//        String errorPattern = ".*Detail: Key \\(core_smiles\\)=\\((.*)\\).*$";
//
//        Pattern p = Pattern.compile(errorPattern);
//        Matcher m = p.matcher(errorMsg2);
//        if(m.matches()){
//            System.out.println(m.group(0));
//            System.out.println(m.group(1));
//        }

//        FrontierDAO dao = FrontierDAO.getInstance();
//        try {
//            Vector<MyCore> allMyCores = dao.getMyCoresBySubstructure("[*]C1=CC=C([*])C=C1");
//            Vector<String> superusers = dao.getSuper_users();
//            for(String superuser:superusers){
//                System.out.println(superuser);
//            }
//            String coreMol = "\n" +
//                    "  Mrv15c1402081710302D          \n" +
//                    "\n" +
//                    "  8  8  0  0  0  0            999 V2000\n" +
//                    "   -1.2385    0.2352    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -1.9530   -0.1773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -1.9530   -1.0023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -1.2385   -1.4148    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -0.5241   -1.0023    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -0.5241   -0.1773    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -1.2385    1.0602    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "   -1.2385   -2.2398    0.0000 R#  0  0  0  0  0  0  0  0  0  0  0  0\n" +
//                    "  1  2  1  0  0  0  0\n" +
//                    "  1  6  2  0  0  0  0\n" +
//                    "  2  3  2  0  0  0  0\n" +
//                    "  3  4  1  0  0  0  0\n" +
//                    "  4  5  2  0  0  0  0\n" +
//                    "  5  6  1  0  0  0  0\n" +
//                    "  1  7  1  0  0  0  0\n" +
//                    "  4  8  1  0  0  0  0\n" +
//                    "M  RGP  2   7   1   8   2\n" +
//                    "M  END\n";
//            String coreSmiles = "[*]C1=CC=C([*])C=C1";
////            System.out.println(dao.getMolFromSmartsFromAldrich("[#6]C(=O)[OD1]",-1));
////            System.out.println(dao.getCASNumberFromStructure("c1ccccc1"));
//            byte[] cdx = Files.toByteArray(new File("/Users/jfeng1/insilico.cdx"));
//            dao.insertMyCoreWrapper("myCore",coreMol,coreSmiles,2,"test",cdx,System.getProperty("user.name"), "My",0);
//        } catch (SQLException e) {
//            System.out.println(e.getMessage());
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

}
