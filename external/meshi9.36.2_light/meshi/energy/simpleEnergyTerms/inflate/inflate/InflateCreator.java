/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.inflate.inflate;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.SegmentList;
import meshi.molecularElements.Protein;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.filters.Filter;
import meshi.util.info.InfoType;

import java.io.File;

public class InflateCreator extends EnergyCreator implements KeyWords {
    private SegmentList segments;
    private InflateType inflateType;
    private File directory = null;
    public final Filter filter;

    public InflateCreator(InflateType inflateType) {
        super(InfoType.INFLATE_ENERGY);
        this.inflateType = inflateType;
        filter = null;
    }

    public InflateCreator(InflateType inflateType, File directory, Filter filter) {
        super(InfoType.INFLATE_ENERGY);
        this.inflateType = inflateType;
        if (inflateType != InflateType.BY_OTHER_MODEL)
            throw new RuntimeException("This constructor my only be used with " + InflateType.BY_OTHER_MODEL);
        this.directory = directory;
        this.filter = filter;
    }


    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        segments = new SegmentList(protein);

        return createEnergyTerm(protein, distanceMatrix, commands, segments);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands, SegmentList segments) {
        Command command = commands.firstWordFilter(infoType).secondWord(RMS_TARGET);
        double rmsTarget = command.thirdWordDouble();
        if (inflateType == InflateType.BY_SEGMENT) {
            EnergyInfoElement info = new EnergyInfoElement(InfoType.INFLATE_BY_SEGMENT, "Inflate by segment", weight);
            term = new InflateBySegments(distanceMatrix, rmsTarget, info, segments);
        } else if (inflateType == InflateType.PER_SEGMENT) {
            EnergyInfoElement info = new EnergyInfoElement(InfoType.INFLATE_PER_SEGMENT, "Inflate per segment", weight);
            term = new InflatePerSegment(distanceMatrix, rmsTarget, info, segments);
        } else if (inflateType == InflateType.BY_OTHER_MODEL) {
            EnergyInfoElement info = new EnergyInfoElement(InfoType.INFLATE_BY_OTHER_MODEL, "Inflate by other model", weight);
            term = new InflateByModel(distanceMatrix, rmsTarget, info, segments, directory, filter);
        } else if (inflateType == InflateType.SIMPLE) {
            EnergyInfoElement info = new EnergyInfoElement(InfoType.INFLATE_ENERGY, "Simple inflate", weight);
            term = new Inflate(distanceMatrix, rmsTarget, info, segments);
        }

        return term;
    }
}