/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry.putH;

import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.angle.*;

import java.util.*;

public class NotEnoughBoundAtomsException extends Exception {
    public NotEnoughBoundAtomsException(String s) {
        super(s);
    }
}
