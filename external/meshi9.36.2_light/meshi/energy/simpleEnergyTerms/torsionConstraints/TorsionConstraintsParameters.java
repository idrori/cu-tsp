package meshi.energy.simpleEnergyTerms.torsionConstraints;

import meshi.energy.Parameters;

/**
 *
 */
public class TorsionConstraintsParameters implements Parameters {
    public final double target;
    public final double std;
    public TorsionConstraintsParameters(double[] targetAndRange) {
        this.target = targetAndRange[0];
        this.std = targetAndRange[1];
    }

}
