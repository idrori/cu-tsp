/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

/**
 * A very liberal filter.
 **/
package meshi.util.filters;

public class AcceptAll implements Filter {
    public boolean accept(Object obj) {
        return true;
    }
}
