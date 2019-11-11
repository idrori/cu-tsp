/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

public class KolDichfin implements Filter {
    public boolean accept(Object obj) {
        if (obj == null) return false;
        return true;
    }
}

