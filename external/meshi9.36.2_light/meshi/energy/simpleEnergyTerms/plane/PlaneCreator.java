/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.plane;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.simpleEnergyTerms.SimpleEnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.info.InfoType;

public class PlaneCreator extends SimpleEnergyCreator implements KeyWords {


    public PlaneCreator() {
        super(InfoType.PLANE);
    }
    public PlaneCreator(CommandList commands) {
        this();
        parametersList = new PlaneParametersList(parametersDirectory(commands) +
                "/" + PLANE_PARAMETERS);
    }

    public AbstractEnergy createEnergyTerm(Protein protein, DistanceMatrix distanceMatrix,
                                           CommandList commands) {
        if (parametersList == null)
            parametersList = new PlaneParametersList(parametersDirectory(commands) +
                    "/" + PLANE_PARAMETERS);

        TorsionList torsionList = TorsionList.createQuickAndDirtyTorsionList(protein, distanceMatrix);
        TorsionList relevantTorsionList = torsionList.filter(new HaveParametersFilter(parametersList));
        // Filtering out a problematic plane torsion in ARG (CD-NE-CZ-NH2)
        TorsionList finalTorsionList = new TorsionList();
        for (Torsion tor : relevantTorsionList)
            if (!(tor.atom1.name().equals("CD") && tor.atom2.name().equals("NE") &&
                    tor.atom3.name().equals("CZ") && tor.atom4.name().equals("NH2")))
                finalTorsionList.add(tor);
        EnergyInfoElement info = new EnergyInfoElement(InfoType.PLANE, "Plane energy", weight);
        term = new PlaneEnergy(finalTorsionList, distanceMatrix, (PlaneParametersList) parametersList, info);
        return term;
    }
}
	    
