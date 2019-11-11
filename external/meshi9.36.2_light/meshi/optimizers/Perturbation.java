/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.optimizers;

import meshi.energy.EvaluationException;
import meshi.util.UpdateableException;


/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 24/02/2010
 * Time: 23:02:24
 * To change this template use File | Settings | File Templates.
 */
public interface Perturbation {
    public void perturb() throws OptimizerException, UpdateableException, EvaluationException;
    public void reset();
    public void restart();
}
