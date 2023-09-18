package com.insilico.application.insilicotools.gui.widget;

/**
 * Created by jfeng1 on 3/10/17.
 */
public class CoreStatus {
    int status_id;
    String status_name;

    public CoreStatus(int status_id, String status_name) {
        this.status_id = status_id;
        this.status_name = status_name;
    }

    public int getStatus_id() {
        return status_id;
    }

    public String getStatus_name() {
        return status_name;
    }

    @Override
    public String toString() {
        return status_name;
    }
}
