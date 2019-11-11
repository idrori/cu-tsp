/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.geometry.GridStatus;

public class UpdateableException extends Exception {
    String comment = "";
    GridStatus gridStatus = null;

    public UpdateableException() {
        super();
    }

    public UpdateableException(GridStatus gridStatus) {
        super(gridStatus.toString() + " " + gridStatus.comment());
        this.gridStatus = gridStatus;
    }
}
