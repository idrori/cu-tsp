/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

public class Equals implements Filter {
    private Object pivot;

    public Equals(Object obj) {
        pivot = obj;
    }

    public boolean accept(Object obj) {
        return pivot.equals(obj);
    }
}
