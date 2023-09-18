package com.insilico.application.insilicotools.data;

import com.google.common.base.Objects;

/**
 * Created by jfeng1 on 8/1/17.
 */
public class Status {
    int status_id;
    String status_name;

    public Status(int status_id, String status_name) {
        this.status_id = status_id;
        this.status_name = status_name;
    }

    public int getStatus_id() {
        return status_id;
    }

    public void setStatus_id(int status_id) {
        this.status_id = status_id;
    }

    public String getStatus_name() {
        return status_name;
    }

    public void setStatus_name(String status_name) {
        this.status_name = status_name;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Status that = (Status) o;
        return status_id == that.status_id &&
                Objects.equal(status_name, that.status_name);
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(status_id, status_name);
    }

    @Override
    public String toString() {
        return status_name;
    }
}
