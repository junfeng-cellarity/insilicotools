package com.insilico.application.insilicotools.data;

import com.google.common.base.Objects;
import openeye.oechem.OEFormat;
import openeye.oechem.OEGraphMol;
import openeye.oechem.oechem;
import openeye.oechem.oemolistream;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;

/**
 * Created by jfeng1 on 9/29/16.
 */
public class SerializableMol implements Serializable{
    int type;
    String molStr;
    String name;
    final static long serialVersionUID = 9078026061736730231L;

    public SerializableMol(String molStr, int type) {
        this.molStr = molStr;
        this.type = type;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public int getType() {
        return type;
    }

    public String getMolStr() {
        return molStr;
    }

    public OEGraphMol getOEMol(){
        OEGraphMol mol = new OEGraphMol();
        oemolistream ifs = new oemolistream();
        ifs.SetFormat(type);
        ifs.openstring(molStr);
        oechem.OEReadMolecule(ifs,mol);
        return mol;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        SerializableMol that = (SerializableMol) o;
        return type == that.type &&
                Objects.equal(molStr, that.molStr);
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(type, molStr);
    }

    public static void main(String[] args) {
        try {
            oechem.OEAddLicenseFromHttp(new URL("http://javelin.corp.My.com:8080/insilico/oe_license.txt"));
        } catch (IOException e) {
            e.printStackTrace();
        }
        ArrayList<SerializableMol> molList = new ArrayList<SerializableMol>();
        SerializableMol mol = new SerializableMol("CCCCCC", OEFormat.SMI);
        SerializableMol mol2 = new SerializableMol("c1ccccc1",OEFormat.SMI);
        molList.add(mol);
        molList.add(mol2);
        try {
//            FileOutputStream fout = new FileOutputStream("/Users/jfeng1/tmp.ser");
//            ObjectOutputStream oop = new ObjectOutputStream(fout);
//            oop.writeObject(molList);
//            oop.close();
//            fout.close();

            FileInputStream fin = new FileInputStream("/Users/jfeng1/tmp.ser");
            ObjectInputStream iip = new ObjectInputStream(fin);
            ArrayList<SerializableMol> molList2 = (ArrayList<SerializableMol>) iip.readObject();
            iip.close();
            fin.close();

            for(SerializableMol m:molList2){
                System.out.println(m.getOEMol().NumAtoms());
            }


        } catch (IOException e) {
            e.printStackTrace();
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        }


    }

}
