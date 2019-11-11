/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.*;

import java.util.*;

public class LineSearchException extends Exception {
    public static final int WEIRD_INPUT_TO_FIND_STEP_LENGTH = 0;
    public static final int NOT_A_DESCENT_DIRECTION = 1;
    public static final int WOLF_CONDITION_NOT_MET = 2;

    public int code;

    public LineSearchException(int code, String msg) {
        super(msg);
        this.code = code;
    }
}



