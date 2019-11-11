/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import meshi.energy.solvation.hydrogenBonds.*;
import meshi.geometry.*;
import meshi.util.Utils;

public class SolvationDistanceAttributePolars extends SolvationDistanceAttribute {
    public double sigmHB;  // The hydrogen bond value (including distance and angular dependence  /*1*/ /*2*/
    public final double saltBridgeFactorA1;     // Salt bridges vs. hydrogen bond strengths /*1*/ /*2*/
    public final double saltBridgeFactorA2;             /*1*//*2*/
    public final int atom1number, atom2number;
    public final AbstractHydrogenBond hydrogenBond;
    public AbstractHydrogenBondList solvationHB;

    public SolvationDistanceAttributePolars(Distance distance, AbstractHydrogenBondList solvateHB,
                                            double saltBridgeFactorA1, double saltBridgeFactorA2, AbstractHydrogenBond hydrogenBond) {
        super(distance, SolvationDistanceType.TWO_POLARS);
        this.solvationHB = solvateHB;
        this.saltBridgeFactorA1 = saltBridgeFactorA1;
        this.saltBridgeFactorA2 = saltBridgeFactorA2;
        this.atom1number = distance.atom1.number;
        this.atom2number = distance.atom2.number;
        this.hydrogenBond = hydrogenBond;
    }

    public void update() {
        // Is this a hydrogen bond or salt bridge?  Possible HB Ahoy!!
        // -----------------------------------------------------------
        AbstractHydrogenBond hb = solvationHB.findBondByPolars(atom1number, atom2number);
        if (hb != null) {
            hb.updateHBvalueAndDerivatives();
            sigmHB = hb.hbVal();
        }
        else
            sigmHB = 0.0;
        return;
    }

    public void debug() {
        AbstractHydrogenBond hb = solvationHB.findBondByPolars(atom1number, atom2number);
        Utils.printDebug(this," "+hb);
    }

}
