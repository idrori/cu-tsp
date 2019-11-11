/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.scoringFunctions;

import meshi.energy.EvaluationException;
import meshi.util.UpdateableException;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfo;

/**

 */
public interface Score  {
    public MeshiInfo score(MeshiInfo energyInfo) ;
}
