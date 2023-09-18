package com.insilico.application.insilicotools.gui.filter.eganEgg;

import com.insilico.application.insilicotools.gui.filter.FilterState;

public class EganEggFilterState extends FilterState {
    private int absState = ABS_95;
    private int bbbState = BBB_90;

    protected final static int ABS_99 = 0;
    protected final static int ABS_95 = 1;
    protected final static int IGNORE = 2;
    protected final static int BBB_99 = 3;
    protected final static int BBB_90 = 4;

    public String isValidState() {
        return absState == IGNORE && bbbState == IGNORE ? "Either ABS or BBB must not be set to 'Ignore'" : null;
    }

    public int getAbsState() {
        return absState;
    }

    public int getBbbState() {
        return bbbState;
    }

    public void setAbsState(final int absState) {
        if(absState != ABS_99 && absState != ABS_95 && absState != IGNORE) {
            throw new IllegalArgumentException("Unknown state value");
        }
        this.absState = absState;
    }

    public void setBbbState(final int bbbState) {
        if(bbbState != BBB_99 && bbbState != BBB_90 && bbbState != IGNORE) {
            throw new IllegalArgumentException("Unknown state value");
        }
        this.bbbState = bbbState;
    }
}