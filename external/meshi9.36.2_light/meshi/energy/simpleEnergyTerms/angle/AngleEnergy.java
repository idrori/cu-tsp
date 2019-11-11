/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.angle;

import meshi.energy.simpleEnergyTerms.*;
import meshi.energy.*;
import meshi.geometry.*;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.MeshiInfoElement;

/**
 * Angle  energy term.
 * Has the general form <b> E = SIGMAi(Ki(Ti-D0i)^2)</b>
 * where <b>Ti</b>
 * is the angle between three consecutive getAtoms, <b> D0i </b>
 * is their expected average angle (depends isOn their types) and <b>Ki</b>
 * is a force constant that again, depends isOn the atom types.<br>
 * This class is used for both calculating the angle-energy term of an energy function
 * and for updating the forces isOn each atom accordingly.<b>
 * It is assumed that the list of angles is constant during the simulation. <b>
 * The numerical method was addapted from Ron Elber's Moil.
 * <p/>
 * Important Note: This energy term has a non-continous point at angle values of 0 or Pi.
 * In order to circumvent these discontinuites we modified the regular parabolic form of this
 * term near the problematic values. At about ~10 degrees (the exact value is hard coded in class
 * AngleEnergyElement) from both 0 and PI the energy starts to climb very steeply so that energetic
 * values of infinity are setResidue to the discontinuous points. We thus hope the simulation could never
 * reach them. On very rare starting condition, however, these problems might never the less be
 * encountered.
 */

public class AngleEnergy extends SimpleEnergyTerm {
    /**
     * The list of angles that needs to be evaluated.
     */
    protected AngleList angleList;
    protected DistanceMatrix distanceMatrix;

    public AngleEnergy(AngleList angleList, DistanceMatrix distanceMatrix,
                       AngleParametersList parametersList, EnergyInfoElement info) {
        super(toArray(distanceMatrix, angleList), parametersList, info);
        this.angleList = angleList;
        this.distanceMatrix = distanceMatrix;
        createElementsList(angleList);
        comment = "Angle";
    }


    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        return new AngleEnergyElement(((Angle) baseElement), parameters, weight);
    }

}
	
    
