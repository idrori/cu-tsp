/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public class MeshiInitException extends Exception {
    private String message;

    public MeshiInitException(String errorMessage) {
        super();
        message = errorMessage + "\n";
    }

    public String message() {
        return message;
    }
}
