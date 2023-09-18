package com.insilico.application.insilicotools.gui.filter;

public class ResultEntry {
    TaskEntry task;
    Throwable err;

    public ResultEntry() {
    }

    public ResultEntry(TaskEntry task) {
        this.task = task;
    }

    public ResultEntry(TaskEntry task, Throwable err) {
        this.task = task;
        this.err = err;
    }

    public final boolean hasError() {
        return err != null;
    }

    public final Throwable getError() {
        return err;
    }

}
