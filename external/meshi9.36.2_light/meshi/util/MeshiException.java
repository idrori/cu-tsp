/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public class MeshiException extends RuntimeException {
    public MeshiException(String errorMessage) {
        super();
        System.err.print("\n" + errorMessage + "\n");
        printStackTrace();
    }
}
