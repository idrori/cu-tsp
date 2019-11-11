package meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranCore;

import meshi.energy.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsionsList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.ResidueTorsionsAttribute;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.parameters.SecondaryStructure;
import meshi.util.MeshiAttribute;

/**
 * Created with IntelliJ IDEA.
 * User: tetyanam
 * Date: 05/02/14
 * Time: 12:46
 * To change this template use File | Settings | File Templates.
 */
public class RamachandranCoreEnergy extends AbstractEnergy implements CompositeTorsionsDefinitions {
    private RamachandranSidechainEnergy ramachandranSidechainEnergy;
    private ResidueTorsionsList residueTorsionsList;

    private static final double minPHI1 = -125*Math.PI/180;
    private static final double maxPHI1 = -30*Math.PI/180;
    private static final double slopePHI1 = 12*Math.PI/180;

    private static final double minPSI1 = -75*Math.PI/180;
    private static final double maxPSI1 = 23*Math.PI/180;
    private static final double slopePSI1 = 12*Math.PI/180;

    private static final double minPHI2 = 40*Math.PI/180;
    private static final double maxPHI2 = 88*Math.PI/180;
    private static final double slopePHI2 = 15*Math.PI/180;

    private static final double minPSI2 = 5*Math.PI/180;
    private static final double maxPSI2 = 75*Math.PI/180;
    private static final double slopePSI2 = 25*Math.PI/180;

    private static final double minPHI3 = -178*Math.PI/180;
    private static final double maxPHI3 = -45*Math.PI/180;
    private static final double slopePHI3 = 12*Math.PI/180;

    private static final double minPSI3 = 90*Math.PI/180;
    private static final double maxPSI3 = 178*Math.PI/180;
    private static final double slopePSI3 = 12*Math.PI/180;

    static BoreCalculator phiBoreForRightHandedHelix;
    static BoreCalculator psiBoreForRightHandedHelix;

    static BoreCalculator phiBoreForLeftHandedHelix;
    static BoreCalculator psiBoreForLeftHandedHelix;

    static BoreCalculator phiBoreForBetaSheet;
    static BoreCalculator psiBoreForBetaSheet;

    double weight;
    private boolean ssFlag;
    //private static final double weight1 = 1, weight2 = 2, weight3 = 1;

    public RamachandranCoreEnergy() {
    }

    public RamachandranCoreEnergy(RamachandranSidechainEnergy ramachandranSidechainEnergy, EnergyInfoElement info, boolean ssFlag, double weight) {
        super(toArray(), info);
        this.ssFlag = ssFlag;
        this.ramachandranSidechainEnergy = ramachandranSidechainEnergy;
        residueTorsionsList = ramachandranSidechainEnergy.residueTorsionsList();
        comment = "ramachandranCore";
        this.weight = weight;
        phiBoreForRightHandedHelix = new BoreCalculator(minPHI1, maxPHI1, slopePHI1);
        psiBoreForRightHandedHelix = new BoreCalculator(minPSI1, maxPSI1, slopePSI1);

        phiBoreForLeftHandedHelix = new BoreCalculator(minPHI2, maxPHI2, slopePHI2);
        psiBoreForLeftHandedHelix = new BoreCalculator(minPSI2, maxPSI2, slopePSI2);

        phiBoreForBetaSheet = new BoreCalculator(minPHI3, maxPHI3, slopePHI3);
        psiBoreForBetaSheet = new BoreCalculator(minPSI3, maxPSI3, slopePSI3);
    }

    public EnergyInfoElement evaluate() {
        return evaluate(false);
    }

    public void evaluateAtoms() {
        evaluate(true);
    }

    public EnergyInfoElement evaluate(boolean evaluateAtoms) {
        double energy;

        Residue residue;
        double phi,psi;
        ResidueTorsionsAttribute rta;

        energy = 0;

        for (ResidueTorsions residueTorsions : residueTorsionsList) {
            if ((residueTorsions.residue().name).equals("PRO") ||
                (residueTorsions.residue().name).equals("GLY")) continue;
            rta = (ResidueTorsionsAttribute) residueTorsions.getAttribute(MeshiAttribute.RESIDUE_TORSIONS_ATTRIBUTE);
            if (rta == null) continue; //if all getAtoms of residue are frozen

            phi = residueTorsions.phi().torsion();
            psi = residueTorsions.psi().torsion();

            if ((!ssFlag) || residueTorsions.residue().getSecondaryStructure() != SecondaryStructure.SHEET)
                energy += makeEnergy(phiBoreForRightHandedHelix, psiBoreForRightHandedHelix, phi, psi, residueTorsions);

            if ((!ssFlag) || residueTorsions.residue().getSecondaryStructure() == SecondaryStructure.COIL);
                energy += makeEnergy(phiBoreForLeftHandedHelix, psiBoreForLeftHandedHelix, phi, psi, residueTorsions);

            if ((!ssFlag) || residueTorsions.residue().getSecondaryStructure() != SecondaryStructure.HELIX)
                energy += makeEnergy(phiBoreForBetaSheet, psiBoreForBetaSheet, phi, psi, residueTorsions);
        }

        if (evaluateAtoms) {
            int nAtoms = 0;
            for (ResidueTorsions residueTorsions : residueTorsionsList) {
                nAtoms += residueTorsions.getResidue().getAtoms().size();
            }
            for (ResidueTorsions residueTorsions : residueTorsionsList) {
                residue = residueTorsions.getResidue();
                AtomList atoms = residue.getAtoms();
                for (Atom atom : atoms)
                    atom.addEnergy(weight*energy / nAtoms);
            }
        }
        energy *= weight;
        info.setValue(energy);
        return info;
    }

    private double makeEnergy(BoreCalculator phiBoreCalculator, BoreCalculator psiBoreCalculator,
                              double phi, double psi, ResidueTorsions residueTorsions){
        double e, phideriv, psideriv;
        double e_phi, de_phi;
        double e_psi, de_psi;

        phiBoreCalculator.calculateBore(phi);
        e_phi = phiBoreCalculator.f();
        de_phi = phiBoreCalculator.df();
        psiBoreCalculator.calculateBore(psi);
        e_psi = psiBoreCalculator.f();
        de_psi = psiBoreCalculator.df();

        e = -e_phi*e_psi;
        phideriv = -e_psi*de_phi;
        psideriv = -e_phi*de_psi;

        residueTorsions.applyForce(PHI, -weight*phideriv);
        residueTorsions.applyForce(PSI, -weight*psideriv);
        return e;
    }

    public void test(TotalEnergy totalEnergy, Atom atom)  {
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
                ramachandranSidechainEnergy.evaluate();
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate().energy();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += EnergyElement.DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
                totalEnergy.update();
                ramachandranSidechainEnergy.evaluate();
            double z2 = evaluate().energy();
            double de = z2 - e1;
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
            if ((z2 == AbstractEnergy.INFINITY) | (z2 == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\nz2 = " + z2);
            if ((analiticalForce == AbstractEnergy.INFINITY) | (analiticalForce == AbstractEnergy.NaN))
                System.out.println("Testing " + this + "\nanaliticalForce = " + analiticalForce);
        }


    }

    }
