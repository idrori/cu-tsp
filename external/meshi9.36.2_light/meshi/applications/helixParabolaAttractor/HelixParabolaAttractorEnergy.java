/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.helixParabolaAttractor;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;
import meshi.energy.TotalEnergy;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionPair;
import meshi.geometry.TorsionPairList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfoElement;

/**
 * Calculates energy value as described in Creator class.
 */
public class HelixParabolaAttractorEnergy extends AbstractEnergy {

    /* find the minimum that attracts our helices */
    private static double phiAttractor = -60 * Math.PI / 180;
    private static double psiAttractor = -40 * Math.PI / 180;
    private static double PI = Math.PI;

    private TorsionPairList helixPairList;


    public HelixParabolaAttractorEnergy(TorsionPairList helixPairList, DistanceMatrix distanceMatrix, EnergyInfoElement info) {
        super(toArray(distanceMatrix, helixPairList), info);
        this.helixPairList = helixPairList;
        this.weight = weight;
        comment = "HelixPara";
    }

    public EnergyInfoElement evaluate() {
        return evaluate(false);
    }

    /*Evaluate attractor energy values. */

    public EnergyInfoElement evaluate(boolean evaluateAtoms) {
        Torsion torsionPhi, torsionPsi; /* the torsion classes for each angle */
        double phi, psi; /* phi, psi values for single torsions pair */
        double phiMinAttractor, psiMinAttractor, phiMinusPI, phiPlusPI, psiMinusPI, psiPlusPI;
        double localEnergy; /* energy for specific pair */
        double derivPhi, derivPsi; /* energy derivative for phi and psi */
        double energy = 0.0; /* total energy achieved so far */
        double energyPerAtom;
        double epsilon = 0.01;
        double alpha = 1;
        double e1phi, e2phi, e3phi, ephi, phiMinAttractor2, phiMinusPI2, phiPlusPI2, phiMinAttractor3, phiMinusPI3, phiPlusPI3, phiMinAttractor4, phiMinusPI4, phiPlusPI4;
        double e1psi, e2psi, e3psi, epsi, psiMinAttractor2, psiMinusPI2, psiPlusPI2, psiMinAttractor3, psiMinusPI3, psiPlusPI3, psiMinAttractor4, psiMinusPI4, psiPlusPI4;

        for (TorsionPair pair : helixPairList) {
            torsionPhi = pair.torsion1();
            torsionPsi = pair.torsion2();

            /* retrieve phi, psi values */
            phi = torsionPhi.torsion();
            psi = torsionPsi.torsion();

            /* calculate local energy and derivation */
            phiMinAttractor = (phi - phiAttractor);
            phiMinusPI = (phi - PI);
            phiPlusPI = (phi + PI);
            phiMinAttractor2 = phiMinAttractor * phiMinAttractor;
            phiMinAttractor3 = phiMinAttractor * phiMinAttractor2;
            phiMinAttractor4 = phiMinAttractor2 * phiMinAttractor2;
            e1phi = (phiMinAttractor4 / (phiMinAttractor4 + alpha)) - 1;
            phiMinusPI2 = phiMinusPI * phiMinusPI;
            phiMinusPI3 = phiMinusPI2 * phiMinusPI;
            phiMinusPI4 = phiMinusPI2 * phiMinusPI2;
            e2phi = phiMinusPI4 / (phiMinusPI4 + epsilon);
            phiPlusPI2 = phiPlusPI * phiPlusPI;
            phiPlusPI3 = phiPlusPI * phiPlusPI2;
            phiPlusPI4 = phiPlusPI2 * phiPlusPI2;
            e3phi = phiPlusPI4 / (phiPlusPI4 + epsilon);

            ephi = e2phi * e1phi * e3phi;


            psiMinAttractor = (psi - psiAttractor);
            psiMinusPI = (psi - PI);
            psiPlusPI = (psi + PI);
            psiMinAttractor2 = psiMinAttractor * psiMinAttractor;
            psiMinAttractor3 = psiMinAttractor * psiMinAttractor2;
            psiMinAttractor4 = psiMinAttractor2 * psiMinAttractor2;
            e1psi = (psiMinAttractor4 / (psiMinAttractor4 + alpha)) - 1;
            psiMinusPI2 = psiMinusPI * psiMinusPI;
            psiMinusPI3 = psiMinusPI * psiMinusPI2;
            psiMinusPI4 = psiMinusPI2 * psiMinusPI2;
            e2psi = psiMinusPI4 / (psiMinusPI4 + epsilon);
            psiPlusPI2 = psiPlusPI * psiPlusPI;
            psiPlusPI3 = psiPlusPI * psiPlusPI2;
            psiPlusPI4 = psiPlusPI2 * psiPlusPI2;
            e3psi = psiPlusPI4 / (psiPlusPI4 + epsilon);

            epsi = e2psi * e1psi * e3psi;

            localEnergy = epsi + ephi + 2;

            derivPhi = (4 * phiMinusPI3 * phiPlusPI3 / ((phiMinusPI4 + epsilon) * (phiPlusPI4 + epsilon))) *
                    (
                            ((alpha * phiMinAttractor3 * phiMinusPI * phiPlusPI) / ((phiMinAttractor4 + alpha) * (phiMinAttractor4 + alpha)))
                                    + (((phiMinAttractor4) / (phiMinAttractor4 + alpha)) - 1) * epsilon * phiPlusPI / (phiMinusPI4 + epsilon)
                                    + (((phiMinAttractor4) / (phiMinAttractor4 + alpha)) - 1) * epsilon * phiMinusPI / (phiPlusPI4 + epsilon)
                    );

            derivPsi = (4 * psiMinusPI3 * psiPlusPI3 / ((psiMinusPI4 + epsilon) * (psiPlusPI4 + epsilon))) *
                    (
                            ((alpha * psiMinAttractor3 * psiMinusPI * psiPlusPI) / ((psiMinAttractor4 + alpha) * (psiMinAttractor4 + alpha)))
                                    + (((psiMinAttractor4) / (psiMinAttractor4 + alpha)) - 1) * epsilon * psiPlusPI / (psiMinusPI4 + epsilon)
                                    + (((psiMinAttractor4) / (psiMinAttractor4 + alpha)) - 1) * epsilon * psiMinusPI / (psiPlusPI4 + epsilon)
                    );


            derivPhi = -1 * derivPhi; /* negation into force */
            derivPsi = -1 * derivPsi; /* negation into force */

            /* apply weight */
            localEnergy = localEnergy * weight;
            energyPerAtom = localEnergy / 8;
            derivPhi = derivPhi * weight;
            derivPsi = derivPsi * weight;

            /* apply to getAtoms */
            if (!torsionPhi.atom1.frozen()) {
                torsionPhi.atom1.addToFx(derivPhi * torsionPhi.dTorsionDx1());
                torsionPhi.atom1.addToFy(derivPhi * torsionPhi.dTorsionDy1());
                torsionPhi.atom1.addToFz(derivPhi * torsionPhi.dTorsionDz1());
            }
            if (!torsionPhi.atom2.frozen()) {
                torsionPhi.atom2.addToFx(derivPhi * torsionPhi.dTorsionDx2());
                torsionPhi.atom2.addToFy(derivPhi * torsionPhi.dTorsionDy2());
                torsionPhi.atom2.addToFz(derivPhi * torsionPhi.dTorsionDz2());
            }
            if (!torsionPhi.atom3.frozen()) {
                torsionPhi.atom3.addToFx(derivPhi * torsionPhi.dTorsionDx3());
                torsionPhi.atom3.addToFy(derivPhi * torsionPhi.dTorsionDy3());
                torsionPhi.atom3.addToFz(derivPhi * torsionPhi.dTorsionDz3());
            }
            if (!torsionPhi.atom4.frozen()) {
                torsionPhi.atom4.addToFx(derivPhi * torsionPhi.dTorsionDx4());
                torsionPhi.atom4.addToFy(derivPhi * torsionPhi.dTorsionDy4());
                torsionPhi.atom4.addToFz(derivPhi * torsionPhi.dTorsionDz4());
            }
            if (!torsionPsi.atom1.frozen()) {
                torsionPsi.atom1.addToFx(derivPsi * torsionPsi.dTorsionDx1());
                torsionPsi.atom1.addToFy(derivPsi * torsionPsi.dTorsionDy1());
                torsionPsi.atom1.addToFz(derivPsi * torsionPsi.dTorsionDz1());
            }
            if (!torsionPsi.atom2.frozen()) {
                torsionPsi.atom2.addToFx(derivPsi * torsionPsi.dTorsionDx2());
                torsionPsi.atom2.addToFy(derivPsi * torsionPsi.dTorsionDy2());
                torsionPsi.atom2.addToFz(derivPsi * torsionPsi.dTorsionDz2());
            }
            if (!torsionPsi.atom3.frozen()) {
                torsionPsi.atom3.addToFx(derivPsi * torsionPsi.dTorsionDx3());
                torsionPsi.atom3.addToFy(derivPsi * torsionPsi.dTorsionDy3());
                torsionPsi.atom3.addToFz(derivPsi * torsionPsi.dTorsionDz3());
            }
            if (!torsionPsi.atom4.frozen()) {
                torsionPsi.atom4.addToFx(derivPsi * torsionPsi.dTorsionDx4());
                torsionPsi.atom4.addToFy(derivPsi * torsionPsi.dTorsionDy4());
                torsionPsi.atom4.addToFz(derivPsi * torsionPsi.dTorsionDz4());
            }

            if (evaluateAtoms) {
                torsionPsi.atom1.addEnergy(energyPerAtom);
                torsionPsi.atom2.addEnergy(energyPerAtom);
                torsionPsi.atom3.addEnergy(energyPerAtom);
                torsionPsi.atom4.addEnergy(energyPerAtom);
                torsionPhi.atom1.addEnergy(energyPerAtom);
                torsionPhi.atom2.addEnergy(energyPerAtom);
                torsionPhi.atom3.addEnergy(energyPerAtom);
                torsionPhi.atom4.addEnergy(energyPerAtom);
            }
            /* update global energy */
            energy += localEnergy;
        }
        info.setValue(energy);
        return info;
    }

    public void evaluateAtoms() {
        evaluate(true);
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        System.out.println("HelixParabolaAttractorEnergy test (not really) fo atom " + atom);
        for (TorsionPair pair : helixPairList) {
            Torsion torsionPhi = pair.torsion1();
            Torsion torsionPsi = pair.torsion2();
            if ((torsionPhi.atom1 == atom) ||
                    (torsionPhi.atom2 == atom) ||
                    (torsionPhi.atom3 == atom) ||
                    (torsionPhi.atom4 == atom) ||
                    (torsionPsi.atom1 == atom) ||
                    (torsionPsi.atom2 == atom) ||
                    (torsionPsi.atom3 == atom) ||
                    (torsionPsi.atom4 == atom))
                System.out.println("Suspected torsion pair:\n" + torsionPhi + "\n" + torsionPsi + "\n----------------\n");
        }
    }


}
