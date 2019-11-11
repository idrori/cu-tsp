/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.simpleEnergyTerms.bond.*;

public interface Classes {
    Class BOND_PARAMETERS_CLASS = (new BondParameters()).getClass();
    Class BOND_ELEMENT_CLASS = (new BondEnergyElement()).getClass();
    //Class LENNARD_JONES_PARAMETERS_CLASS = (new LennardJonesParameters()).getClass();    
    //Class LENNARD_JONES_ELEMENT_CLASS    = (new LennardJonesEnergyElement()).getClass();
}
