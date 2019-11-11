/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.energy.*;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

public class RgEnergy extends AbstractEnergy {
    private RgInfoElement rgInfo;
    public Atom testAtom = null;
    private final RgParameters n_RGhSS = RgParameters.n_RGhSS;
    private final RgParameters RGhSS_RGhCoilCentered = RgParameters.RGhSS_RGhCoilCentered;
    private final RgParameters RGhSS_RGbSScentered = RgParameters.RGhSS_RGbSScentered;
    private final RgParameters RGhSS_RGbCoilCentered = RgParameters.RGhSS_RGbCoilCentered;
    private final RgParameters RGhSS_RGcSScentered = RgParameters.RGhSS_RGcSScentered;
    private final RgParameters RGhSS_RGcCoilCentered = RgParameters.RGhSS_RGcCoilCentered;

    private final double logN;
    private final LogRadiusOfGyration logNonPolarSsRG;
    private final LogRadiusOfGyration logNonPolarCoilRG;
    private final LogRadiusOfGyration logBackboneSsRG;
    private final LogRadiusOfGyration logBackboneCoilRG;
    private final LogRadiusOfGyration logPolarSsRG;
    private final LogRadiusOfGyration logPolarCoilRG;
    private final CenterOfMass centerOfMass, backboneCenterOfMass;
    private final RadiusOfGyration nonPolarRG, backboneRG, polarRG;


    public RgEnergy(RgInfoElement info) {
        this(getDummyUpdatables(), new LogRadiusOfGyration[6],0,info);
        off();
    }
    public RgEnergy(Updateable[] updateables, LogRadiusOfGyration[] logRgArray, double logN, RgInfoElement info) {
        super(updateables,new EnergyInfoElement(InfoType.RG,InfoType.RG.tag), EnergyType.NON_DIFFERENTIAL);
        comment = "RG";
        this.logN = logN;
        logNonPolarSsRG = logRgArray[0];
        logNonPolarCoilRG = logRgArray[1];
        logBackboneSsRG = logRgArray[2];
        logBackboneCoilRG = logRgArray[3];
        logPolarSsRG = logRgArray[4];
        logPolarCoilRG = logRgArray[5];
        if (updateables[0] instanceof CenterOfMass)
            centerOfMass = (CenterOfMass) updateables[0];
        else
            centerOfMass = null;
        if (updateables[1] instanceof CenterOfMass)
            backboneCenterOfMass = (CenterOfMass) updateables[1];
        else backboneCenterOfMass = null;
        if (logNonPolarSsRG != null)
            nonPolarRG = logNonPolarSsRG.rg;
        else nonPolarRG = null;
        if (logBackboneSsRG != null)
            backboneRG = logBackboneSsRG.rg;
        else backboneRG = null;
        if (logPolarSsRG != null)
            polarRG = logPolarSsRG.rg;
        else polarRG = null;
        setInfoElements();
        rgInfo = (RgInfoElement) info;
        weight = rgInfo.all;
    }

    public void evaluateAtoms() {
    }

    public EnergyInfoElement evaluate() {
        if (logNonPolarSsRG == null) off();
        if (!isOn()) return info;

        double energy = 0;

        rgInfo.N_RG_HSS.setValue(evaluateN_RGhSS(rgInfo.RGlogaritmicWeight));
        rgInfo.E_HSS_HCOIL.setValue(evaluateRGratio(centerOfMass, logNonPolarSsRG, logNonPolarCoilRG, RGhSS_RGhCoilCentered, testAtom, rgInfo.RGratioWeight, rgInfo.HSS_HCOIL));
        rgInfo.E_HSS_BSS.setValue(evaluateRGratio(centerOfMass, logNonPolarSsRG, logBackboneSsRG, RGhSS_RGbSScentered, testAtom, rgInfo.RGratioWeight, rgInfo.HSS_BSS));
        rgInfo.E_HSS_BCOIL.setValue(evaluateRGratio(centerOfMass, logNonPolarSsRG, logBackboneCoilRG, RGhSS_RGbCoilCentered, testAtom, rgInfo.RGratioWeight, rgInfo.HSS_BCOIL));
        rgInfo.E_HSS_PSS.setValue(evaluateRGratio(centerOfMass, logNonPolarSsRG, logPolarSsRG, RGhSS_RGcSScentered, testAtom, rgInfo.RGratioWeight, rgInfo.HSS_PSS));
        rgInfo.E_HSS_PCOIL.setValue(evaluateRGratio(centerOfMass, logNonPolarSsRG, logPolarCoilRG, RGhSS_RGcCoilCentered, testAtom, rgInfo.RGratioWeight, rgInfo.HSS_PCOIL));
        rgInfo.HSS.setValue(evaluateLinearRG(centerOfMass, nonPolarRG, rgInfo.RGnonPolarWeight));
        rgInfo.BSS.setValue(evaluateLinearRG(backboneCenterOfMass, backboneRG, rgInfo.RGbackboneWeight));
        rgInfo.PSS.setValue(evaluateLinearRG(centerOfMass, polarRG, rgInfo.RGpolarWeight));
        energy += ((Double) rgInfo.N_RG_HSS.getValue()).doubleValue();
        energy += ((Double) rgInfo.E_HSS_HCOIL.getValue()).doubleValue();
        energy += ((Double) rgInfo.E_HSS_BSS.getValue()).doubleValue();
        energy += ((Double) rgInfo.E_HSS_BCOIL.getValue()).doubleValue();
        energy += ((Double) rgInfo.E_HSS_PSS.getValue()).doubleValue();
        energy += ((Double) rgInfo.E_HSS_PCOIL.getValue()).doubleValue();
        energy += ((Double) rgInfo.HSS.getValue()).doubleValue();
        energy += ((Double) rgInfo.BSS.getValue()).doubleValue();
        energy += ((Double) rgInfo.PSS.getValue()).doubleValue();
        info.getChildren().clear();
        info.getChildren().add(rgInfo.N_RG_HSS);
        info.getChildren().add(rgInfo.E_HSS_HCOIL);
        info.getChildren().add(rgInfo.E_HSS_BSS);
        info.getChildren().add(rgInfo.E_HSS_BCOIL);
        info.getChildren().add(rgInfo.E_HSS_PSS);
        info.getChildren().add(rgInfo.E_HSS_PCOIL);
        info.getChildren().add(rgInfo.HSS);
        info.getChildren().add(rgInfo.BSS);
        info.getChildren().add(rgInfo.PSS);
        info.getChildren().add(rgInfo.HSS_HCOIL);
        info.getChildren().add(rgInfo.HSS_BSS);
        info.getChildren().add(rgInfo.HSS_BCOIL);
        info.getChildren().add(rgInfo.HSS_PSS);
        info.getChildren().add(rgInfo.HSS_PCOIL);
        info.setValue(energy);
        return info;
    }


    private double evaluateLinearRG(CenterOfMass centerOfMass, RadiusOfGyration radiusOfGyration, double factor) {
        double rg = radiusOfGyration.radiusOfGyration();
        double energy = weight * factor * rg;
        Atom atom;
        int size = radiusOfGyration.size();
        int sizeCM = centerOfMass.size();

        for (int i = 0; i < size; i++) {
            atom = radiusOfGyration.atomAt(i);
            atom.addToFx(-weight * factor * radiusOfGyration.dRGdX(i));
            atom.addToFy(-weight * factor * radiusOfGyration.dRGdY(i));
            atom.addToFz(-weight * factor * radiusOfGyration.dRGdZ(i));
        }
        for (int i = 0; i < sizeCM; i++) {
            atom = centerOfMass.atomAt(i);
            atom.addToFx(-weight * factor * radiusOfGyration.dRGdCMX());
            atom.addToFy(-weight * factor * radiusOfGyration.dRGdCMy());
            atom.addToFz(-weight * factor * radiusOfGyration.dRGdCMz());
        }
        return energy;


    }

    private double evaluateN_RGhSS(double myWeight) {
        double energy = 0;

        double logRG = 0;
        if (logNonPolarSsRG != null)
            logRG = logNonPolarSsRG.logRG();

        double dRatio = 1 / (n_RGhSS.alpha * logN + 1);
        double ratio = (logRG - n_RGhSS.beta) * dRatio;

        double diffU = ratio - n_RGhSS.upperLimit;
        double diffL = ratio - n_RGhSS.lowerLimit;
        double diff;
        if ((diffU <= 0) && (diffL >= 0)) return 0;
        if (diffU > 0) diff = diffU;
        else diff = diffL;
        diff = 100 * diff;
        double diff3 = diff * diff * diff;
        double diff4 = diff3 * diff;
        energy = weight * diff4 * myWeight;

        int size = logNonPolarSsRG.size();
        Atom atom;
        for (int i = 0; i < size; i++) {
            atom = logNonPolarSsRG.atomAt(i);
            atom.addToFx(-weight * myWeight * 400 * diff3 * dRatio * logNonPolarSsRG.dLogRGdX(i));
            atom.addToFy(-weight * myWeight * 400 * diff3 * dRatio * logNonPolarSsRG.dLogRGdY(i));
            atom.addToFz(-weight * myWeight * 400 * diff3 * dRatio * logNonPolarSsRG.dLogRGdZ(i));
        }
        return energy;
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        System.out.println(" testing RgEnergy using " + atom);
        if (!on) {
            System.out.println("" + this + " is off");
            return;
        }
        testAtom = atom;
        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for (int i = 0; i < 3; i++) {
                totalEnergy.update();
            double x = coordinates[i][0];
            coordinates[i][1] = 0;
            double e1 = evaluate().energy();
            double analiticalForce = coordinates[i][1];
            coordinates[i][0] += EnergyElement.DX;
            // Whatever should be updated ( such as distance matrix torsion list etc. )
                totalEnergy.update();
            double e2 = evaluate().energy();
            double de = e2 - e1;
            double numericalForce = -de / EnergyElement.DX;
            coordinates[i][0] -= EnergyElement.DX;
                totalEnergy.update();

            double diff = Math.abs(analiticalForce - numericalForce);

            if (Math.abs(analiticalForce - numericalForce) > 0.00001) {
                //if ((2*diff/(Math.abs(analiticalForce)+Math.abs(numericalForce)+VERY_SMALL)) > relativeDiffTolerance){
                System.out.println("Testing " + this);
                System.out.println("Atom[" + atom.number() + "]." + EnergyElement.XYZ.charAt(i) + " = " + x);
                System.out.println("Analytical force = " + analiticalForce);
                System.out.println("Numerical force  = " + numericalForce);


                System.out.println("diff = " + diff + "\n" +
                        "tolerance = 2*diff/(|analiticalForce| + |numericalForce|+VERY_SMALL) = " +
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

    private double evaluateRGratio(CenterOfMass centerOfMass, LogRadiusOfGyration logRG1, LogRadiusOfGyration logRG2, RgParameters parameters,
                                   Atom testAtom, double myWeight, MeshiInfo infoElement) {
        double energy = 0;
        int size1 = logRG1.size();
        int size2 = logRG2.size();
        int sizeCM = centerOfMass.size();
        double logRg1 = logRG1.logRG();


        double logRg2 = logRG2.logRG();

        double ratio = (logRg2 - parameters.beta) / (logRg1 * parameters.alpha + 1);
        infoElement.setValue(ratio);
        double diffU = ratio - parameters.upperLimit;
        double diffL = ratio - parameters.lowerLimit;
        double diff;
        if ((diffU <= 0) && (diffL >= 0)) return 0;
        if (diffU > 0) diff = diffU;
        else diff = diffL;
        double diff3 = diff * diff * diff;
        double diff4 = diff3 * diff;
        energy = weight * myWeight * diff4;


        //
        double dEnergy = weight * 4 * diff3;

        double dEnergyDrg1 = -dEnergy * (logRG2.logRG() - parameters.beta) / (logRG1.logRG() * logRG1.logRG() * parameters.alpha);
        double dEnergyDrg2 = dEnergy / (logRG1.logRG() * parameters.alpha);
        double dEnergyDcmX = dEnergyDrg1 * logRG1.dLogRGdCMx() + dEnergyDrg2 * logRG2.dLogRGdCMx();
        double dEnergyDcmY = dEnergyDrg1 * logRG1.dLogRGdCMy() + dEnergyDrg2 * logRG2.dLogRGdCMy();
        double dEnergyDcmZ = dEnergyDrg1 * logRG1.dLogRGdCMz() + dEnergyDrg2 * logRG2.dLogRGdCMz();

        Atom atom;
        for (int i = 0; i < size1; i++) {
            atom = logRG1.atomAt(i);
            atom.addToFx(-dEnergyDrg1 * logRG1.dLogRGdX(i));
            atom.addToFy(-dEnergyDrg1 * logRG1.dLogRGdY(i));
            atom.addToFz(-dEnergyDrg1 * logRG1.dLogRGdZ(i));
        }
        for (int i = 0; i < size2; i++) {
            atom = logRG2.atomAt(i);
            atom.addToFx(-dEnergyDrg2 * logRG2.dLogRGdX(i));
            atom.addToFy(-dEnergyDrg2 * logRG2.dLogRGdY(i));
            atom.addToFz(-dEnergyDrg2 * logRG2.dLogRGdZ(i));
        }
        for (int i = 0; i < sizeCM; i++) {
            atom = centerOfMass.atomAt(i);
            atom.addToFx(-dEnergyDcmX);
            atom.addToFy(-dEnergyDcmY);
            atom.addToFz(-dEnergyDcmZ);
        }
        return energy;
    }


public void scaleWeight(double factor) {
    super.scaleWeight(factor);
}
    public void setInfoElements() {

    }

    public void test(Atom atom) {
        testLinearRG(centerOfMass, polarRG, rgInfo.RGpolarWeight,atom);
    }

    private void testLinearRG(CenterOfMass centerOfMass, RadiusOfGyration radiusOfGyration, double factor, Atom atom) {
              System.out.println(" testLinearRG(CenterOfMass centerOfMass, RadiusOfGyration radiusOfGyration, double factor, Atom atom) ");
            double delta = 0.00000001;
            atom.setFx(0);
            atom.setFy(0);
            atom.setFz(0);
            double rg = radiusOfGyration.radiusOfGyration();
            double energy = weight * factor * rg;
            System.out.println("Analytical fx="+atom.fx()+"  fy="+atom.fy()+"  fz="+atom.fz());
            atom.addToX(delta);
            rg = radiusOfGyration.radiusOfGyration();
            double  energy1 = weight * factor * rg;
            System.out.println("Numerical X = "+(energy-energy1)/delta);
            atom.addToX(-delta);
            atom.addToY(delta);
            rg = radiusOfGyration.radiusOfGyration();
             energy1 = weight * factor * rg;
             System.out.println("Numerical Y = "+(energy-energy1)/delta);
             atom.addToY(-delta);
             atom.addToZ(delta);
             rg = radiusOfGyration.radiusOfGyration();
             energy1 = weight * factor * rg;
             System.out.println("Numerical Z = "+(energy-energy1)/delta);
             atom.addToZ(-delta);

             // radiusOfGyration.test(atom);
             // centerOfMass.test(atom);
        }

    private static Updateable[] getDummyUpdatables() {
        Updateable[] out = new Updateable[2];
        out[0] = new DummyUpdatable();
        out[1] = new DummyUpdatable();
        return out;
    }
    private static class DummyUpdatable implements Updateable {
        public void update(int numberOfUpdates) throws UpdateableException {}

        public void setNumberOfUpdates(int numberOfUpdates){}

    }

}
      