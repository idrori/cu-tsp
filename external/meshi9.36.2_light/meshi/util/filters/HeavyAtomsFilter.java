/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.filters;

import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

public class HeavyAtomsFilter implements Filter {
    public boolean accept(Object obj) {
        Atom atom = (Atom) obj;
        if (atom.nowhere()) throw new RuntimeException("Does not know  what to do with a \"nowhere\" atom " + atom);
        if (atom.type().isHydrogen()) return false;
        return true;
    }
}
