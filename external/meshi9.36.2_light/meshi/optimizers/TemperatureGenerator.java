/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;
/*import meshi.util.*;
import meshi.molecularElements.*;
import meshi.energy.*;
import meshi.energy.inflate.*;
import meshi.energy.distanceConstraints.*;
*/

import java.util.*;

public class TemperatureGenerator {
    public final double initialTemperature;
    public final double finalTemperature;
    public final double de;
    private double currentTemperature;

    public TemperatureGenerator(double initialTemperature, double finalTemperature, int maxSteps) {
        this.initialTemperature = initialTemperature;
        this.finalTemperature = finalTemperature;
        de = (finalTemperature - initialTemperature) / maxSteps;
        currentTemperature = initialTemperature;
    }

    public double currentTemperature() {
        return currentTemperature;
    }

    public double next() {
        currentTemperature += de;
        return currentTemperature;
    }
}
