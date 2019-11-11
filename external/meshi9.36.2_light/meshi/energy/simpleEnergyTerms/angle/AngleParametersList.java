/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.angle;

import meshi.util.*;
import meshi.energy.simpleEnergyTerms.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.util.filters.Filter;
import meshi.util.file.*;
import meshi.energy.*;

import java.util.*;

public class AngleParametersList extends ParametersList {
    public AngleParametersList(String parametersFileName) {
        super(parametersFileName, true);
    }

    public Parameters createParameters(String line) {
        return new AngleParameters(line);
    }

    public Parameters parameters(Object obj) {
        Angle angle = (Angle) obj;
        Parameters key = new AngleParameters(angle.atom1.type(), angle.atom2.type(), angle.atom3.type());
        return getParameters(key);
    }
}
