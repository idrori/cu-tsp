/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

/**
 *
 */
public enum GridStatus {
    OK, FAILED_AXIS_CHANGED_TOO_FAST, FAILED_GRID_TOO_LARGE;

    private String comment;

    private GridStatus() {
        comment = "No comment";
    }

    public String comment() {
        return comment;
    }

    public GridStatus setComment(String comment) {
        this.comment = comment;
        return this;
    }
}
