/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.geometry.*;
import meshi.parameters.*;
import meshi.energy.solvate.hydrogenBonds.*;

public class SolvateDistanceAttributeCreator {
    public static SolvateDistanceAttribute create(Distance dis, SolvateParametersList parameters, AbstractHydrogenBondList solvateHB) {
        // Does this distance involve a hydrogen ??
        // ----------------------------------------
        if (dis.atom1().type().isHydrogen() || dis.atom2().type().isHydrogen())
            return SolvateEnergy.NOT_PATRICIPATING_DISTANCE;

        // Does this distance involves two polar atoms ??
        // ----------------------------------------------

        if ((dis.atom1().type().isOxygen() || dis.atom1().type().isNitrogen() || (dis.atom1().type() == AtomType.CSG)) &&
                (dis.atom2().type().isOxygen() || dis.atom2().type().isNitrogen() || (dis.atom2().type() == AtomType.CSG))) {


            //---------------------
            double saltBridgeFactorA1, saltBridgeFactorA2;
            if ((((dis.atom1().type() == AtomType.KNZ) || (dis.atom1().type() == AtomType.RNH) || (dis.atom1().type() == AtomType.TRN)) &&
                    ((dis.atom2().type() == AtomType.DOD) || (dis.atom2().type() == AtomType.EOE) || (dis.atom2().type() == AtomType.TRO))) ||
                    (((dis.atom2().type() == AtomType.KNZ) || (dis.atom2().type() == AtomType.RNH) || (dis.atom2().type() == AtomType.TRN)) &&
                            ((dis.atom1().type() == AtomType.DOD) || (dis.atom1().type() == AtomType.EOE) || (dis.atom1().type() == AtomType.TRO)))) {  // This is a salt bridge
                if (dis.atom1().type() == AtomType.DOD)
                    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ASP_OD;
                else if (dis.atom1().type() == AtomType.EOE)
                    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_GLU_OE;
                else if (dis.atom1().type() == AtomType.TRO)
                    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRO;
                else if (dis.atom1().type() == AtomType.KNZ)
                    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_LYS_NZ;
                else if (dis.atom1().type() == AtomType.RNH)
                    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ARG_NH;
                else if (dis.atom1().type() == AtomType.TRN)
                    saltBridgeFactorA1 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRN;
                else saltBridgeFactorA1 = 1.0;
                if (dis.atom2().type() == AtomType.DOD)
                    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ASP_OD;
                else if (dis.atom2().type() == AtomType.EOE)
                    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_GLU_OE;
                else if (dis.atom2().type() == AtomType.TRO)
                    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRO;
                else if (dis.atom2().type() == AtomType.KNZ)
                    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_LYS_NZ;
                else if (dis.atom2().type() == AtomType.RNH)
                    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_ARG_NH;
                else if (dis.atom2().type() == AtomType.TRN)
                    saltBridgeFactorA2 = SolvateEnergy.SALT_BRIDGE_STRENGTH_TRN;
                else saltBridgeFactorA2 = 1.0;
            } else {
                saltBridgeFactorA1 = 1.0;
                saltBridgeFactorA2 = 1.0;
            }
            return new SolvateDistanceAttributeBetweenPolars(dis, solvateHB,
                    saltBridgeFactorA1, saltBridgeFactorA2);
        }
        //------------------------------------
        int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        //=================================

        SigmaParameters sigmaParameters12 = parameters.sigmaParameters[TsaiAtomicType1][TsaiAtomicType2];
        SigmaParameters sigmaParameters21 = parameters.sigmaParameters[TsaiAtomicType2][TsaiAtomicType1];

        if ((dis.atom1.type().isCarbon() || (dis.atom1.type() == AtomType.MSD)) &&
                (dis.atom1.type().isCarbon() || (dis.atom1.type() == AtomType.MSD))) {
            return new SolvateDistanceAttributeWithNonPolar(dis, sigmaParameters12, sigmaParameters21);
        }
        if (dis.atom1.type().isCarbon() || (dis.atom1.type() == AtomType.MSD))
            return new SolvateDistanceAttributeNonPolarPolar(dis, sigmaParameters21);

        return new SolvateDistanceAttributePolarNonPolar(dis, sigmaParameters12);
    }
}