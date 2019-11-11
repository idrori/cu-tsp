/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import meshi.geometry.*;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.*;
import meshi.energy.solvation.hydrogenBonds.*;

public class SolvationDistanceAttributeCreator {
    public SolvationDistanceAttribute create(Distance dis, SolvateParametersList parameters, AbstractHydrogenBondList solvateHB) {
        // Does this distance involve a hydrogen ??
        // ----------------------------------------
        if (dis.atom1().type().isHydrogen() || dis.atom2().type().isHydrogen())
            return SolvationEnergy.NOT_PATRICIPATING_DISTANCE;
        if ((dis.atom1().type() == AtomType.TRO)|| (dis.atom1().type() == AtomType.TRN)||(dis.atom1().type() == AtomType.TRC) ||
                (dis.atom2().type() == AtomType.TRO)|| (dis.atom2().type() == AtomType.TRN)||(dis.atom2().type() == AtomType.TRC))
            return SolvationEnergy.NOT_PATRICIPATING_DISTANCE;
        int seqSeparation = dis.atom1().residueNumber()-dis.atom2().residueNumber();
        if ((seqSeparation <= 1) && (seqSeparation >= -1)) return SolvationEnergy.NOT_PATRICIPATING_DISTANCE;

        if (!acceptCoil()) {
            if ((dis.atom1().residue().getSecondaryStructure().equals(SecondaryStructure.COIL)) ||
                    (dis.atom2().residue().getSecondaryStructure().equals(SecondaryStructure.COIL)))
                return SolvationEnergy.NOT_PATRICIPATING_DISTANCE;
        }
        // Does this distance involves two polar getAtoms ??
        // ----------------------------------------------

        if ((dis.atom1().type().isOxygen() || dis.atom1().type().isNitrogen() || (dis.atom1().type() == AtomType.CSG)) &&
                (dis.atom2().type().isOxygen() || dis.atom2().type().isNitrogen() || (dis.atom2().type() == AtomType.CSG)) &&
                (dis.atom1().type() != AtomType.PN) && (dis.atom2().type() !=  AtomType.PN)) {

            //---------------------
            double saltBridgeFactorA1, saltBridgeFactorA2;
            if ((((dis.atom1().type() == AtomType.KNZ) || (dis.atom1().type() == AtomType.RNH) ) &&
                    ((dis.atom2().type() == AtomType.DOD) || (dis.atom2().type() == AtomType.EOE) )) ||
                    (((dis.atom2().type() == AtomType.KNZ) || (dis.atom2().type() == AtomType.RNH) ) &&
                            ((dis.atom1().type() == AtomType.DOD) || (dis.atom1().type() == AtomType.EOE)))) {  // This is a salt bridge
                if (dis.atom1().type() == AtomType.DOD)
                    saltBridgeFactorA1 = SolvationEnergy.SALT_BRIDGE_STRENGTH_ASP_OD;
                else if (dis.atom1().type() == AtomType.EOE)
                    saltBridgeFactorA1 = SolvationEnergy.SALT_BRIDGE_STRENGTH_GLU_OE;
                else if (dis.atom1().type() == AtomType.TRO)
                    saltBridgeFactorA1 = SolvationEnergy.SALT_BRIDGE_STRENGTH_TRO;
                else if (dis.atom1().type() == AtomType.KNZ)
                    saltBridgeFactorA1 = SolvationEnergy.SALT_BRIDGE_STRENGTH_LYS_NZ;
                else if (dis.atom1().type() == AtomType.RNH)
                    saltBridgeFactorA1 = SolvationEnergy.SALT_BRIDGE_STRENGTH_ARG_NH;
                else if (dis.atom1().type() == AtomType.TRN)
                    saltBridgeFactorA1 = SolvationEnergy.SALT_BRIDGE_STRENGTH_TRN;
                else saltBridgeFactorA1 = 1.0;
                if (dis.atom2().type() == AtomType.DOD)
                    saltBridgeFactorA2 = SolvationEnergy.SALT_BRIDGE_STRENGTH_ASP_OD;
                else if (dis.atom2().type() == AtomType.EOE)
                    saltBridgeFactorA2 = SolvationEnergy.SALT_BRIDGE_STRENGTH_GLU_OE;
                else if (dis.atom2().type() == AtomType.TRO)
                    saltBridgeFactorA2 = SolvationEnergy.SALT_BRIDGE_STRENGTH_TRO;
                else if (dis.atom2().type() == AtomType.KNZ)
                    saltBridgeFactorA2 = SolvationEnergy.SALT_BRIDGE_STRENGTH_LYS_NZ;
                else if (dis.atom2().type() == AtomType.RNH)
                    saltBridgeFactorA2 = SolvationEnergy.SALT_BRIDGE_STRENGTH_ARG_NH;
                else if (dis.atom2().type() == AtomType.TRN)
                    saltBridgeFactorA2 = SolvationEnergy.SALT_BRIDGE_STRENGTH_TRN;
                else saltBridgeFactorA2 = 1.0;
            } else {
                saltBridgeFactorA1 = 1.0;
                saltBridgeFactorA2 = 1.0;
            }
            AbstractHydrogenBond hydrogenBond =  solvateHB.findBondByPolars(dis.atom1(),dis.atom2());
            if (hydrogenBond == null) return SolvationEnergy.NOT_PATRICIPATING_DISTANCE;
            return new SolvationDistanceAttributePolars(dis, solvateHB,
                    saltBridgeFactorA1, saltBridgeFactorA2,hydrogenBond);
        }
        //------------------------------------
        int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        //=================================

        SigmaParameters sigmaParameters12 = parameters.sigmaParameters[TsaiAtomicType1][TsaiAtomicType2];
        SigmaParameters sigmaParameters21 = parameters.sigmaParameters[TsaiAtomicType2][TsaiAtomicType1];

        SolvationDistanceType type;
        if (dis.atom1.type().isCarbon() || (dis.atom1.type() == AtomType.MSD) || (dis.atom1.type() == AtomType.PN)) {
            if (dis.atom2.type().isCarbon() || (dis.atom2.type() == AtomType.MSD) || (dis.atom2.type() == AtomType.PN))
                type = SolvationDistanceType.TWO_HYDROPHOBIC;
            else
                type = SolvationDistanceType.HYDROPHOBIC_POLAR;
        }
        else type = SolvationDistanceType.POLAR_HYDROPHOBIC;

        return new SolvationDistanceAttribute(dis, type,sigmaParameters12, sigmaParameters21);
    }

    public boolean acceptCoil() {return true;}
    public boolean isDisulfiedBonded(Atom atom) {
        if (!(atom.type() == AtomType.CSG)) return false;
        for (Atom bonded : atom.bonded()) {
            if (bonded.type() == AtomType.CSG) return true;
        }
        return false;
    }
}