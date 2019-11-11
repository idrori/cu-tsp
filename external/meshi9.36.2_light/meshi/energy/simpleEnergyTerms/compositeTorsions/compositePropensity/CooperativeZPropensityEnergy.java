/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.compositeTorsions.compositePropensity;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.ResidueType;
import meshi.util.MeshiAttribute;

/**
 * Created by IntelliJ IDEA.
 * User: tetyanam
 * Date: 09/06/2009
 * Time: 10:09:03
 * To change this template use File | Settings | File Templates.
 */
public class CooperativeZPropensityEnergy extends AbstractEnergy implements CompositeTorsionsDefinitions {
    private static final double[] mean = new double[ResidueType.values().length];
    private static final double[] std = new double[ResidueType.values().length];
    private static final double WIDTH = 10;
    private static final double HEIHT = 1;
    private static double M = (HEIHT + 1) * WIDTH;

    private CompositePropensityEnergy propensityEnergy;
    private ResidueTorsionsList residueTorsionsList;
    int nResidues;

    public CooperativeZPropensityEnergy() {
    }

    public CooperativeZPropensityEnergy(CompositePropensityEnergy propensityEnergy, EnergyInfoElement info, CooperativeZPropensityParameters parameters) {
        super(toArray(), info);
        this.propensityEnergy = propensityEnergy;
        residueTorsionsList = propensityEnergy.residueTorsionsList;
        comment = "CoopZPropensity";
        for (ResidueType type : ResidueType.values()) {
            mean[type.ordinal()] = parameters.mean[type.ordinal()];
            std[type.ordinal()] = parameters.std[type.ordinal()];
        }
    }

    public EnergyInfoElement evaluate() {
        return evaluate(false);
    }

    public void evaluateAtoms() {
        evaluate(true);
    }

    public EnergyInfoElement evaluate(boolean evaluateAtoms) {

        double energy, z, z3, z2, zMinusEpsilon, zMinusEpsilon2;
        double derivative, factor;
        Residue residue;
        ResidueType type;
        z = 0;
        nResidues = residueTorsionsList.size();
        for (ResidueTorsions rt : residueTorsionsList) {
            type = rt.getResidueType();
            if (std[type.ordinal()] != 0)
                z += (rt.energy() - mean[type.ordinal()]) / std[type.ordinal()];
        }
        z /= nResidues;

        z2 = z * z;
        energy = weight * (M / WIDTH - (M + z2) / (z2 + WIDTH));
        derivative = (-1) * weight * 2 * z * (WIDTH - M) / (z2 + WIDTH) / (z2 + WIDTH) / nResidues;

        if (evaluateAtoms) {
            int nAtoms = 0;
            for (ResidueTorsions residueTorsions : residueTorsionsList) {
                nAtoms += residueTorsions.getResidue().getAtoms().size();
            }
            for (ResidueTorsions residueTorsions : residueTorsionsList) {
                residue = residueTorsions.getResidue();
                AtomList atoms = residue.getAtoms();
                for (Atom atom : atoms)
                    atom.addEnergy((weight * energy) / nAtoms);
            }
        }

        ResidueTorsionsPropensityAttribute rta;
        for (ResidueTorsions residueTorsions : residueTorsionsList) {
            rta = (ResidueTorsionsPropensityAttribute) residueTorsions.getAttribute(MeshiAttribute.RESIDUE_TORSIONS_ATTRIBUTE);
            if (rta == null) continue; //if all getAtoms of residue are frozen
            residue = residueTorsions.getResidue();
            type = residue.type();
            if (std[type.ordinal()] == 0)
                throw new RuntimeException("Something is weird in Parameters or residueTorsionsList of" + this);
            factor = derivative / std[type.ordinal()];
            if (factor != 0) {
                residueTorsions.applyForce(PHI, -rta.phi_deriv * factor);
                residueTorsions.applyForce(PSI, -rta.psi_deriv * factor);
            }
        }
        info.setValue(energy);        //*40
        return info;
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (!on) {
            System.out.println("" + this + " is off");
            return;
        }
        System.out.println("Testing " + this + " using " + atom);
        if (atom == null)
            throw new RuntimeException("Cannot test " + this);

        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for (int i = 0; i < 3; i++) {
                totalEnergy.update();
            propensityEnergy.evaluate();
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate().energy();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += EnergyElement.DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
                totalEnergy.update();
            propensityEnergy.evaluate();
            double e2 = evaluate().energy();
            double de = e2 - e1;
            double numericalForce = -de / EnergyElement.DX;
            coordinates[i][0] -= EnergyElement.DX;
                totalEnergy.update();

            double diff = Math.abs(analiticalForce - numericalForce);

            if ((2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + EnergyElement.VERY_SMALL)) > EnergyElement.relativeDiffTolerance) {
                System.out.println("Testing " + this);
                System.out.println("Atom[" + atom.number() + "]." + EnergyElement.XYZ.charAt(i) + " = " + x);
                System.out.println("Analytical force = " + analiticalForce);
                System.out.println("Numerical force  = " + numericalForce);

                System.out.println("diff = " + diff + "\n" +
                        "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+EnergyElement.VERY_SMALL) = " +
                        2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + EnergyElement.VERY_SMALL));
                System.out.println();
            }
            if ((e1 == AbstractEnergy.INFINITY) | (e1 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\ne1 = " + e1);
            if ((e2 == AbstractEnergy.INFINITY) | (e2 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\ne2 = " + e2);
            if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\nanaliticalForce = " + analiticalForce);
        }

    }


}
