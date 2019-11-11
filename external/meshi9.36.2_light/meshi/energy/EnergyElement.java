/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.UpdateableException;
import meshi.util.formats.Fdouble;

/**
 * The essence of an energy term - its application to a minimal setResidue of getAtoms. The energy term is a summation of
 * (typically) many such elements.
 */
public abstract class EnergyElement {
    protected boolean frozen;
    protected AtomList atoms;
    // protected double[][] coordinates = new double[3][];
    public final static double DX = 1e-7;  /* should be roughly the sqrare root of the machine's relative precision.
					In java this is about sqrt(1e-15) */
    // should be 3-4 orders of magnitude higer than DX.
    public final static double relativeDiffTolerance = 0.0005;
    public final static double VERY_SMALL = Math.exp(-15);
    public final static String XYZ = "xyz";
    protected static Fdouble dFormatStd = Fdouble.STANDARD;
    protected static Fdouble dFormatSrt = Fdouble.SHORT;

    public EnergyElement() {
        atoms = null;
        frozen = false;
    }

    public final boolean frozen() {
        return frozen;
    }

    public final AtomList atoms() {
        return atoms;
    }

    public boolean updateFrozen() {
        frozen = true;
        for (Atom atom : atoms)
            frozen = atom.frozen() & frozen;
        return frozen;
    }

    public abstract double evaluate();

    public double evaluateAtoms() {
        double e = evaluate();
        int nAtoms = atoms.size();
        for (Atom atom : atoms)
            atom.addEnergy(e / nAtoms);
        return e;
    }

    protected abstract void setAtoms();

    public void test(TotalEnergy totalEnergy, Atom atom) {
        if (atoms == null)
            throw new RuntimeException("Cannot test " + this + "\n" + "No getAtoms defined");
        if (atoms.whereIs(atom) < 0)
            return;

        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for (int i = 0; i < 3; i++) {
            totalEnergy.update();

            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
                totalEnergy.update();
            double e2 = evaluate();
            double de = e2 - e1;
            double numericalForce = -de / DX;
            coordinates[i][0] -= DX;
            totalEnergy.update();

            double diff = Math.abs(analiticalForce - numericalForce);

            if (Math.abs(analiticalForce - numericalForce) > 0.00001) {
                //if ((2*diff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+VERY_SMALL)) > relativeDiffTolerance){
                System.out.println("Testing " + this);
                System.out.println("Atom[" + atom.number() + "]." + XYZ.charAt(i) + " = " + x);
                System.out.println("Analytical force = " + analiticalForce);
                System.out.println("Numerical force  = " + numericalForce);
                System.out.println("Weird atoms are:");
                for (Atom a : atoms)
                    System.out.println(a);

                System.out.println("diff = " + diff + "\n" +
                        "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+VERY_SMALL) = " +
                        2 * diff / (Math.abs(analiticalForce) + Math.abs(numericalForce) + VERY_SMALL));
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

