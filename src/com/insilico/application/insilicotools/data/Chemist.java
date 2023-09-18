package com.insilico.application.insilicotools.data;

public class Chemist {
    int chemist_id;
    int project_id;
    String chemist_name;

    public Chemist(int chemist_id, int project_id, String chemist_name) {
        this.chemist_id = chemist_id;
        this.project_id = project_id;
        this.chemist_name = chemist_name;
    }

    public int getChemist_id() {
        return chemist_id;
    }

    public void setChemist_id(int chemist_id) {
        this.chemist_id = chemist_id;
    }

    public int getProject_id() {
        return project_id;
    }

    public void setProject_id(int project_id) {
        this.project_id = project_id;
    }

    public String getChemist_name() {
        return chemist_name;
    }

    public void setChemist_name(String chemist_name) {
        this.chemist_name = chemist_name;
    }

    @Override
    public String toString() {
        return chemist_name;
    }
}
