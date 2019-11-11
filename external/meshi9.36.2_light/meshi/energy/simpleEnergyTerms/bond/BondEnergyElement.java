/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.simpleEnergyTerms.bond;

import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.energy.*;

public class BondEnergyElement extends EnergyElement {
    public static final double TOLERANCE_FACTOR = 4;
    public static final double MAX_ENERGY = 1000;
    protected Atom atom1, atom2;
    protected int number1, number2;
    protected Distance distance;
    protected double target, force, force2, limit;
    private boolean on;
    double weight;
    private double tolerance;
    private boolean simpleFlag;

    public BondEnergyElement() {
    }

    public BondEnergyElement(AtomPair atomPair, Parameters parameters,
                             DistanceMatrix distanceMatrix, double weight,boolean simpleFlag){
        this.weight = weight;
        this.simpleFlag = simpleFlag;
        atom1 = atomPair.atom1();
        atom2 = atomPair.atom2();
        setAtoms();
        int atom1Number = atomPair.atom1Number();
        int atom2Number = atomPair.atom2Number();
        distance = distanceMatrix.distance(atom1Number, atom2Number);
        try {
            if (distance == null) distance = new DistanceMirror(distanceMatrix.distance(atom2Number, atom1Number));
        }
        catch (RuntimeException ex) {
            System.out.println("A problem while building a BondEnergyElement\n" +
                    "atom1 = " + atom1 + "\n" +
                    "atom2 = " + atom2 + "\n" +
                    "distanceMatrix.distance(atom1Number, atom2Number) = " + distanceMatrix.distance(atom1Number, atom2Number) + "\n" +
                    "distanceMatrix.distance(atom2Number, atom1Number) = " + distanceMatrix.distance(atom2Number, atom1Number) + "\n" +
                    "distanceMatrix.bondedList().size() = " + distanceMatrix.bondedList().size());

            for (DistanceList distanceRow : distanceMatrix.bondedList())
                for (Distance distance : distanceRow)
                System.out.println("bonded = " + distance);
            throw ex;
        }

        target = ((BondParameters) parameters).target;
        force = ((BondParameters) parameters).force * weight;
        force2 = ((BondParameters) parameters).force2 * weight;
        updateFrozen();
        on = true;
        tolerance = target*target*TOLERANCE_FACTOR;
        double initD = distance.distance();
        if (initD > tolerance) limit = initD;
        else limit = tolerance;
        if (simpleFlag) limit = 1000;
    }

    public void scaleWeight(double factor) {
        weight *= factor;
    }

/*    public void test(TotalEnergy totalEnergy, Atom atom) throws EvaluationException{
        System.out.println("Testing "+this+"\n"+
                "Distance = "+distance);
        super.test(totalEnergy,atom);
    }  */

    public boolean updateFrozen() {
        super.updateFrozen();
        frozen = frozen || (distance == null);
        return frozen;
    }

    protected void setAtoms() {
        atoms = new AtomList(atom1.molecularSystem);
        atoms.add(atom1);
        atoms.add(atom2);
    }

    public double evaluate() {//throws BondException{
        double d, d2, d3;
        double deDd;
        double deDx;
        double deDy;
        double deDz;
        double energy;
        double e; // The energy when the distance is normal
        double ePlusMax;
        double eTimesMax;
        double ALPHA = 0.0001;
        double d2PlusAlpha;
        double d2PlusAlpha2;
        double dis = distance.distance();
        if (frozen()) return 0.0;
        if (!on) return 0.0;
        d = dis - target;
        d2 = d * d;

//        if (d2 > tolerance) throw new BondException("The bond \n"+distance.atom1()+
//                                                    "\n"+distance.atom2()+"\nis "+dis+" - too long for stable simulation");
        e = d2 * force;
        if (simpleFlag) {
            energy = e;
            deDd = 2*d*force;
        }
        else {
            ePlusMax = e+MAX_ENERGY;
            eTimesMax = e*MAX_ENERGY;
            energy = eTimesMax/ePlusMax;
            deDd = (MAX_ENERGY*d * force2/ePlusMax)*(1-e/ePlusMax);
        }

        if (dis > limit) {
            energy += (dis-limit)*(dis-limit)*1000;
            deDd   +=  2000*(dis-limit);
        }
        if (dis < 0.2) {
            energy += (dis - 0.2) * (dis - 0.2) * 100000;
            deDd += 200000 * (dis - 0.2);
        }

        deDx = deDd * distance.dDistanceDx();
        deDy = deDd * distance.dDistanceDy();
        deDz = deDd * distance.dDistanceDz();
        if (!atom1.frozen()) {
            atom1.addToFx(-1 * deDx); // force = -derivative
            atom1.addToFy(-1 * deDy);
            atom1.addToFz(-1 * deDz);
        }
        if (!atom2.frozen()) {
            atom2.addToFx(deDx);
            atom2.addToFy(deDy);
            atom2.addToFz(deDz);
        }
        if ((!(energy> 0)) && (!(energy == 0)) &&(!(energy<0)) )
            throw new RuntimeException(" weird bond energy element "+energy+" "+dis+" "+limit+" "+MAX_ENERGY+" "+force+"\n"+atom1+"\n"+atom2);
        return energy;
    }

    public String toString() {
        double energy;
//        try {
            energy = evaluate();
//        } catch (EvaluationException ex){
//            System.out.println(ex);
//            energy = 99999999;
//        }
        return ("BondEnergyElement target = " + dFormatSrt.f(target) + " force = " + dFormatSrt.f(force) + " distance = " +
                dFormatSrt.f(distance.distance()) + " energy = " + dFormatSrt.f(energy) + "\n" +
                atom1.verbose(1) + "\n" + atom2.verbose(1));
    }

    public Atom atom1() {
        return atom1;
    }

    public Atom atom2() {
        return atom2;
    }

    public boolean isOn() {
        return on;
    }

    public void turnOn() {
        on = true;
    }

    public void turnOff() {
        on = false;
    }

}
