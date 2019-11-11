/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;

public class CbFilter implements Filter {
    public static final CbFilter filter = new CbFilter();

    public boolean accept(Object obj) {
        if (((Atom) obj).name().equals("CB")) return true;
        return false;
    }
}
