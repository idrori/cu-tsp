/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;

public class ChargedAtomFilter implements Filter {
    public static final ChargedAtomFilter filter = new ChargedAtomFilter();

    public boolean accept(Object obj) {
        return ((Atom) obj).type().isCharged();
    }
}