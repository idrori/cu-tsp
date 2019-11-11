/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology;

import meshi.util.filters.Filter;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;

public class CaAtomFilter implements Filter {
    public boolean accept(Object obj) {
        return ((Atom) obj).name.equals("CA");
    }
}
