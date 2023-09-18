package com.insilico.application.insilicotools.gui.lims;

import java.util.Date;

/**
 * Created by jfeng1 on 3/23/17.
 */
public class Batch {
    int batch_id;
    String chemist;
    String batch_name;
    Date date;

    public Batch(int batch_id, String chemist, String batch_name, Date date) {
        this.batch_id = batch_id;
        this.chemist = chemist;
        this.batch_name = batch_name;
        this.date = date;
    }

    public Date getDate() {
        return date;
    }

    public int getBatch_id() {
        return batch_id;
    }

    public void setBatch_id(int batch_id) {
        this.batch_id = batch_id;
    }

    public String getChemist() {
        return chemist;
    }

    public void setChemist(String chemist) {
        this.chemist = chemist;
    }

    public String getBatchName() {
        return batch_name;
    }

    public void setBatchName(String batch_name) {
        this.batch_name = batch_name;
    }

    @Override
    public int hashCode() {
        return batch_id;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Batch that = (Batch) o;
        return batch_id == that.batch_id;
    }

    @Override
    public String toString() {
        return "Batch{" +
                "batch_id=" + batch_id +
                ", chemist='" + chemist + '\'' +
                ", batch_name='" + batch_name + '\'' +
                '}';
    }
}
