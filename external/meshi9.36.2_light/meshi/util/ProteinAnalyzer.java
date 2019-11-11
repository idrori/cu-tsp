/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.EvaluationException;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentException;
import meshi.util.info.ProteinInfo;
import meshi.util.info.ProteinInfoList;

/**
 */
public interface ProteinAnalyzer {
    public ProteinInfo analyze(Protein protein) throws UpdateableException, EvaluationException,AlignmentException;
}
