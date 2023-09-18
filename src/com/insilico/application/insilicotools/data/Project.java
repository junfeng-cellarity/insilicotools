package com.insilico.application.insilicotools.data;

import com.google.common.base.Objects;

/**
 * Created by jfeng1 on 8/1/17.
 */
public class Project {
    int project_id;
    String project_name;

    public Project(int project_id, String project_name) {
        this.project_id = project_id;
        this.project_name = project_name;
    }

    public int getProject_id() {
        return project_id;
    }

    public void setProject_id(int project_id) {
        this.project_id = project_id;
    }

    public String getProject_name() {
        return project_name;
    }

    public void setProject_name(String project_name) {
        this.project_name = project_name;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Project that = (Project) o;
        return project_id == that.project_id &&
                Objects.equal(project_name, that.project_name);
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(project_id, project_name);
    }

    @Override
    public String toString() {
        return project_name;
    }
}
