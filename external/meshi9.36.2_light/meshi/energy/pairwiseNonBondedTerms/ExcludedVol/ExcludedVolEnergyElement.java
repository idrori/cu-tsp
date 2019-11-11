/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.ExcludedVol;

import meshi.energy.pairwiseNonBondedTerms.NonBondedEnergyElement;
import meshi.geometry.Distance;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;


public class ExcludedVolEnergyElement extends NonBondedEnergyElement {
    protected DistanceMatrix distanceMatrix;
    protected Atom atom1, atom2;
    protected Distance distance;
    int type;
    protected double sigma, C1, C2, Rfac;
    protected boolean frozen;
    protected double dEdD;
    protected double dEdX;
    protected double dEdY;
    protected double dEdZ;
    protected double energy;
    protected double rMax;
    protected ExcludedVolParametersList parametersList;



    public ExcludedVolEnergyElement(ExcludedVolParametersList parametersList,
                                    DistanceMatrix distanceMatrix,
                                    int type,
                                    double weight,
                                    double Rfac) {
        this.parametersList = parametersList;
        this.type = type;
        this.weight = weight;
        this.distanceMatrix = distanceMatrix;
        this.rMax = distanceMatrix.rMax();
        this.Rfac = Rfac;
        atoms = new AtomList(distanceMatrix.molecularSystem);
        //----------------------
    }

    protected void setAtoms() {
        throw new RuntimeException("setAtoms() may not be used by ExcludedVolEnergyElement for " +
                "efficiency.");
    }


    public void set(Object obj) {
        distance = (Distance) obj;

        atom1 = distance.atom1();
        atom2 = distance.atom2();
        atoms.clear();
        atoms.add(atom1);
        atoms.add(atom2);


        ExcludedVolParameters parameters = (ExcludedVolParameters) parametersList.parameters(distance);
        sigma = parameters.sigma;
        C1 = parameters.C1;
        C2 = parameters.C2;
        if (sigma > rMax)
            throw new RuntimeException("Excluded Volume: sigma=" + sigma + " and it is larger " +
                    "than rMax=" + rMax);
    }


    public double evaluate() {
        double dEdD;
        double DMinusSig;
        double D3;
        double localSigma = sigma;
        if (Rfac < 0.99) {
            if (atom1.residueNumber() == atom2.residueNumber())
                localSigma = 0;
            else
                localSigma = sigma * Rfac;
        }
        double dis = distance.distance();
        if (dis > localSigma) {
            energy = 0;
        } else {
            DMinusSig = (dis - localSigma);
            if (type == 0) {
                if ((atom1.name().length() == 1) || (atom2.name().length() == 1) || atom1.name().equals("CA")
                        || atom2.name().equals("CA") || atom1.name().equals("CB") || atom2.name().equals("CB")) {
                    D3 = C2 * DMinusSig * DMinusSig * DMinusSig * weight;
                    energy = D3 * DMinusSig;
                    dEdD = 4 * D3;
                } else {
                    D3 = C1 * DMinusSig * weight;
                    energy = D3 * DMinusSig;
                    dEdD = 2 * D3;
                }
                if (!atom1.frozen()) {
                    atom1.addToFx(-1 * dEdD * distance.dDistanceDx()); // force = -derivative
                    atom1.addToFy(-1 * dEdD * distance.dDistanceDy());
                    atom1.addToFz(-1 * dEdD * distance.dDistanceDz());
                }
                if (!atom2.frozen()) {
                    atom2.addToFx(dEdD * distance.dDistanceDx());
                    atom2.addToFy(dEdD * distance.dDistanceDy());
                    atom2.addToFz(dEdD * distance.dDistanceDz());
                }
            } else if (type == 1) {
                if (((atom1.name().length() == 1) || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
                        (atom2.name().equals("CB") || atom2.name().equals("CA") || (atom2.name().length() == 1))) {
                    D3 = 1.0 / (0.1 * 0.1) * DMinusSig * weight;
                    energy = D3 * DMinusSig;
                    dEdD = 2 * D3;
                } else {
                    energy = 0.0;
                    dEdD = 0.0;
                }
                if (!atom1.frozen()) {
                    atom1.addToFx(-1 * dEdD * distance.dDistanceDx()); // force = -derivative
                    atom1.addToFy(-1 * dEdD * distance.dDistanceDy());
                    atom1.addToFz(-1 * dEdD * distance.dDistanceDz());
                }
                if (!atom2.frozen()) {
                    atom2.addToFx(dEdD * distance.dDistanceDx());
                    atom2.addToFy(dEdD * distance.dDistanceDy());
                    atom2.addToFz(dEdD * distance.dDistanceDz());
                }
            } else if (type == 2) {
                if (((atom1.name().length() == 1) || atom1.name().equals("CA") || atom1.name().equals("CB")) &&
                        (atom2.name().equals("CB") || atom2.name().equals("CA") || (atom2.name().length() == 1))) {
                    D3 = DMinusSig * weight;
                    energy = D3 * DMinusSig;
                    dEdD = 2 * D3;
                } else {
                    energy = 0.0;
                    dEdD = 0.0;
                }
                if (!atom1.frozen()) {
                    atom1.addToFx(-1 * dEdD * distance.dDistanceDx()); // force = -derivative
                    atom1.addToFy(-1 * dEdD * distance.dDistanceDy());
                    atom1.addToFz(-1 * dEdD * distance.dDistanceDz());
                }
                if (!atom2.frozen()) {
                    atom2.addToFx(dEdD * distance.dDistanceDx());
                    atom2.addToFy(dEdD * distance.dDistanceDy());
                    atom2.addToFz(dEdD * distance.dDistanceDz());
                }
            } else
                throw new RuntimeException("An unknown excluded volume type");
        }
        return energy;
    }


    public String toString() {
        if ((atom1 == null) & (atom2 == null)) return "ExcludedVolumeEnergyElement - atoms not setResidue yet";
        if ((atom1 == null) | (atom2 == null)) throw new RuntimeException("This is weird\n" +
                "atom1 = " + atom1 + "\n" +
                "atom2 = " + atom2);
        double dis = distance.distance();
        return ("ExcludedVolumeEnergyElement sigma = " + sigma + " Distance = " +
                dFormatSrt.f(dis) + " rMax = " + rMax + "\n" + atom1.verbose(1) + "\n" + atom2.verbose(1));
    }
}
