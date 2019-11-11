package meshi.features;

import meshi.energy.EnergyInfoElement;
import meshi.energy.EvaluationException;

/**
 * Created by chen on 02/12/2015.
 */
public interface Feature {
    public abstract FeatureInfoElement evaluate() throws EvaluationException;

}
