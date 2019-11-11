package meshi.energy.simpleEnergyTerms.DistanceConstraints;

import meshi.energy.Parameters;

/**
 *
 */
public class DistanceConstraintsParameters implements Parameters {
    public final double target;
    public final double range;
    public DistanceConstraintsParameters(double[] targetAndRange) {
        this.target = targetAndRange[0];
        this.range = targetAndRange[1];
    }

    public DistanceConstraintsParameters(double target) {
        this.target = target;
        this.range = 0;
    }


}
