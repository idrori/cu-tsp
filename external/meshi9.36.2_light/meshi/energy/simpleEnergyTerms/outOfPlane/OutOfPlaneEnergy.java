/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.outOfPlane;

import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfoElement;

/**
 * OutOfPlane energy term.
 */
public class OutOfPlaneEnergy extends SimpleEnergyTerm {
    /**
     * The constructor associates any outOfPlane with its parameters.
     */
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public OutOfPlaneEnergy() {
    }

    public OutOfPlaneEnergy(DistanceMatrix distanceMatrix,
                            TorsionList torsionList,
                            OutOfPlaneParametersList parametersList,
                            EnergyInfoElement info) {
        super(toArray(distanceMatrix, torsionList), parametersList, info);
        this.torsionList = torsionList;
        this.distanceMatrix = distanceMatrix;
        createElementsList(torsionList);
        comment = "OOPlane";

    }


    public Parameters createParameters(String line) {
        return new OutOfPlaneParameters(line);
    }


    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        return new OutOfPlaneEnergyElement(((Torsion) baseElement), parameters, weight);
    }

    public Parameters getKey(Object baseElement) {
        Torsion torsion = (Torsion) baseElement;
        return new OutOfPlaneParameters(torsion.atom1.type(), torsion.atom2.type(),
                torsion.atom3.type(), torsion.atom4.type());
    }

    public void handleMissingParameters(Object obj) {
    }
}	
	
