/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

public class BackboneFilter implements Filter {
    public static final BackboneFilter filter = new BackboneFilter();

    public boolean accept(Object obj) {
        return ((Atom) obj).type().backbone();
    }
}
