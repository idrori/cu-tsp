/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

import java.util.*;

import meshi.util.filters.*;
import meshi.util.*;

public class GridCell extends ArrayList<AtomCore> {
    public GridCell(int capacity) {
        super(capacity);
    }
}
