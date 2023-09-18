package com.insilico.application.insilicotools.gui.filter.dnd;

import java.util.concurrent.atomic.AtomicBoolean;

public class DragAndDropLock {
    private static AtomicBoolean locked = new AtomicBoolean(false);
    private static AtomicBoolean startedDnD = new AtomicBoolean(false);

    public static boolean isLocked() {
        return locked.get();
    }

    public static void setLocked(boolean isLocked) {
        locked.set(isLocked);
    }

    public static boolean isDragAndDropStarted() {
        return startedDnD.get();
    }

    public static void setDragAndDropStarted(boolean isLocked) {
        startedDnD.set(isLocked);
    }
}