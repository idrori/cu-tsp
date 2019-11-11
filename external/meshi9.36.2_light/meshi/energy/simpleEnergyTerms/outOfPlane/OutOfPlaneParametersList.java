/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.outOfPlane;

import meshi.energy.simpleEnergyTerms.*;
import meshi.util.*;
import meshi.molecularElements.*;
import meshi.parameters.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;

import java.util.*;

public class OutOfPlaneParametersList extends ParametersList {
    public OutOfPlaneParametersList(String parametersFileName) {
        super(parametersFileName, true);
    }

    public Parameters createParameters(String line) {
        return new OutOfPlaneParameters(line);
    }

    public Parameters parameters(Object obj) {
        Torsion torsion = (Torsion) obj;
        Parameters key = new OutOfPlaneParameters(torsion.atom1.type(), torsion.atom2.type(),
                torsion.atom3.type(), torsion.atom4.type());
        return getParameters(key);
    }
}
