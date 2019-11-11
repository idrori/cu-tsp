/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.*;
import meshi.util.filters.*;

public class GoodBonds implements Filter {
    public boolean accept(Object obj) {
        AtomPair pair = (AtomPair) obj;
        if ((!pair.atom1().active()) | (!pair.atom2().active())) return false;
        return true;
    }
}
