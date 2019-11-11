/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.bond;

import meshi.util.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.*;

import java.util.*;

/**
 * ParametersList for BondEnergy. The list is sortable.
 */
public class BondParametersList extends ParametersList {
    public BondParametersList(String parametersFileName) {
        super(parametersFileName, true);
    }

    public Parameters createParameters(String line) {
        return new BondParameters(line);
    }

    /**
     * Returns BondParameters object (which specify target-distance and force-constant)
     * for a pair of getAtoms).
     */
    public Parameters parameters(Object obj) {
        AtomPair pair = (AtomPair) obj;
        Parameters key = new BondParameters(pair.largeType(), pair.smallType());
        return getParameters(key);
    }

    public String toString() {
        return " BondParametersList ";
    }
}