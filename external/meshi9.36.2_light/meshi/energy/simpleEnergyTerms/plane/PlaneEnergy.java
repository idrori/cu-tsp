/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.plane;

import meshi.energy.EnergyElement;
import meshi.energy.EnergyInfoElement;
import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.SimpleEnergyTerm;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionList;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;

/**
 * Plane energy term.
 */
public class PlaneEnergy extends SimpleEnergyTerm {
    protected TorsionList torsionList;
    protected DistanceMatrix distanceMatrix;

    public PlaneEnergy() {
    }

    public PlaneEnergy(TorsionList torsionList, DistanceMatrix distanceMatrix,
                       PlaneParametersList parametersList, EnergyInfoElement info) {
        super(toArray(distanceMatrix, torsionList), parametersList, info);
        this.torsionList = torsionList;
        this.distanceMatrix = distanceMatrix;
        createElementsList(torsionList);
        comment = "Plane";

    }

    public void scaleWeight(double factor) {
        for (EnergyElement element:elementsList)
            ((PlaneEnergyElement)element).scaleWeight(factor);
    }


    public EnergyElement createElement(Object baseElement, Parameters parameters) {
        PlaneEnergyElement element;
        element =  new PlaneEnergyElement(((Torsion) baseElement), parameters, weight);
        if ((element.isomer == PlaneParameters.Isomer.CIS) & (element.torsion.cosTorsion()< 0)) {
                throw new PlaneException(element.torsion,element.isomer);
        }
        if ((element.isomer == PlaneParameters.Isomer.TRANS) & (element.torsion.cosTorsion()>= 0)) {
            if (NQspecialCase(element) ||
                backboneElement(element))   element.off();// This is the result of a very bad initial structure after minimization with bond and angle it will anyway get closer to the expected torsion angle.
                else throw new PlaneException(element.torsion,element.isomer);

        }
        return element;
    }

    private static boolean NQspecialCase(PlaneEnergyElement element){
        if ((element.atom1.type() != AtomType.NCB &&
             element.atom1.type() != AtomType.QCG )) return false;
        if ((element.atom2.type() != AtomType.NND &&
             element.atom2.type() != AtomType.QNE )) return false;
        if ((element.atom3.type() != AtomType.NCG &&
             element.atom3.type() != AtomType.QCD )) return false;
        if ((element.atom4.type() != AtomType.NOD &&
             element.atom4.type() != AtomType.QOE )) return false;
        return true;
    }

    private boolean backboneElement(PlaneEnergyElement element) {
        for (Atom atom : element.atoms())
           if ((!atom.isBackbone()) &&  (atom.type()!= AtomType.PCD)) return false;
        return true;
    }
}
