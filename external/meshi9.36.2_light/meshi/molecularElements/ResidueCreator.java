/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

public interface ResidueCreator {
    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem);
}
