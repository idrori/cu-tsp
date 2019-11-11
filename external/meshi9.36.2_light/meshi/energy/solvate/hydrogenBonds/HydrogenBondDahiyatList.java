/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate.hydrogenBonds;

import meshi.energy.solvate.SolvateParametersList;
import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;

/**
 * This class describe a list of hydrogen bonds in a protein. The hydrogen bonds are described as in:
 * Dahiyat et al. Protein Sci. 1996 May;5(5):895-903. See the "HydrogenBondDahiyat" class for a much more
 * detailed description.
 * <p/>
 * The construcor requires a parameter that points to an instance of "SolvateParametersList". From this instance
 * the distance dependence of the hydrogen bonds are extracted. See the heading documentation of the
 * SolvateParametersList class for more details.
 * <p/>
 * <p/>
 * This parameter 30 determines the refresh rate of the HB vector. Once every this number of updates, broken hydrogen
 * bonds that are no longer in the non-bonded-list are removed from the vector. This number should be >>1 so that
 * the refresh does not impend the updates too much, but 30 was just a thumb figure.
 */

public class HydrogenBondDahiyatList extends AbstractHydrogenBondList {

    /**
     * This array stores pointers to the base atoms of every non-hydrogen polar atom in the protein.
     * The indexation is through the atom number (field)
     */
    protected Atom[] baseAtom;
    /**
     * Since currently the hydrogen bond list is only used by the solvate energy, we use the parameters
     * for the hydrogen bond from this class. Later, this can be changed.
     */
    SolvateParametersList parameters;


    public HydrogenBondDahiyatList() {
        throw new RuntimeException("\nERROR: without parameters the hydrogen bonds cannot be formed.\n");
    }

    public HydrogenBondDahiyatList(DistanceMatrix dm, AtomList atomList, SolvateParametersList parameters) {
        super(dm, atomList, 30 /* DEFAULT_REFRESH_VECTOR_RATE */);
        this.parameters = parameters;
        try {
            update(false, 0);
        }
        catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("\nAn error occur while creating the Dahiyat hydrogen bond list.\n\n" + e + "\n\n");
        }
    }


    /**
     * This method updates the pointers in the baseAtom array. For non-polar atoms they should remain null.
     * For any polar atom (O,N or S in Cys) we calculate the base atom that participate in the definition of
     * the hydrogen bond angle. This base atom is the attached hydrogen (if present), or the heavy atom to which
     * the polar atom is attached (when the hydrogen is not present).
     */
    protected void buildSpecificStructures() {
        baseAtom = new Atom[maxAtomNum + 1];
        Atom atom1, atom2;
        for (int c1 = 0; c1 < atomList.size(); c1++) {
            atom1 = atomList.atomAt(c1);
            for (int c2 = 0; c2 < atom1.bonded().size(); c2++) {
                atom2 = atom1.bonded().atomAt(c2);
                // Treating OXYGEN atoms.
                // This is easy because the oxygens in proteins
                // are always tied to one atom only.
                if (atom1.type().isOxygen())
                    baseAtom[lut[atom1.number()]] = atom2;
                // Treating NITROGEN atoms.
                // First, the case of amides in glutamines and asparagines
                if (!atom2.type().isHydrogen() && ((atom1.type() == AtomType.NND) || (atom1.type() == AtomType.QNE)))
                    baseAtom[lut[atom1.number()]] = atom2;
                // Second, the case of amides without explicit H attached
                if ((atom1.type() == AtomType.KNZ) || (atom1.type() == AtomType.RNH) || (atom1.type() == AtomType.TRN))
                    baseAtom[lut[atom1.number()]] = atom2;
                // Third , regular H attached
                if (atom2.type().isHydrogen() && ((atom1.type() == AtomType.HND) || (atom1.type() == AtomType.HNE) ||
                        (atom1.type() == AtomType.RNE) || (atom1.type() == AtomType.WNE) ||
                        (atom1.type() == AtomType.AN) ||
                        (atom1.type() == AtomType.CN) ||
                        (atom1.type() == AtomType.DN) ||
                        (atom1.type() == AtomType.EN) ||
                        (atom1.type() == AtomType.FN) ||
                        (atom1.type() == AtomType.GN) ||
                        (atom1.type() == AtomType.HN) ||
                        (atom1.type() == AtomType.IN) ||
                        (atom1.type() == AtomType.KN) ||
                        (atom1.type() == AtomType.LN) ||
                        (atom1.type() == AtomType.MN) ||
                        (atom1.type() == AtomType.NN) ||
                        (atom1.type() == AtomType.QN) ||
                        (atom1.type() == AtomType.RN) ||
                        (atom1.type() == AtomType.SN) ||
                        (atom1.type() == AtomType.TN) ||
                        (atom1.type() == AtomType.VN) ||
                        (atom1.type() == AtomType.WN) ||
                        (atom1.type() == AtomType.YN)))
                    baseAtom[lut[atom1.number()]] = atom2;
                // Treating the SULFUR atoms of Cystines.
                if (atom1.type() == AtomType.CSG)
                    baseAtom[lut[atom1.number()]] = atom2;
            }
        }
    }

    /**
     * Creating a Dahiyat-like hydrogen-bond.
     */
    protected AbstractHydrogenBond createHBfromPolars(Atom atom1, Atom atom2) {
        if ((atom1.type() == AtomType.MSD) || (atom2.type() == AtomType.MSD) || (atom1.type() == AtomType.PN) || (atom2.type() == AtomType.PN) ||
                atom1.nowhere() || atom2.nowhere() || (baseAtom[lut[atom1.number()]]).nowhere() || (baseAtom[lut[atom2.number()]]).nowhere())
            return null;

        AtomList tmpList = new AtomList(atom1.molecularSystem);
        tmpList.add(atom1);
        tmpList.add(atom2);
        tmpList.add(baseAtom[lut[atom1.number()]]);
        tmpList.add(baseAtom[lut[atom2.number()]]);

        int TsaiAtomicType1 = parameters.atomicTypeConverter[atom1.type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        int TsaiAtomicType2 = parameters.atomicTypeConverter[atom2.type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        double end = parameters.hbParameters[TsaiAtomicType1][TsaiAtomicType2].end;
        double p1 = parameters.hbParameters[TsaiAtomicType1][TsaiAtomicType2].p1;
        double p2 = parameters.hbParameters[TsaiAtomicType1][TsaiAtomicType2].p2;
        double valAtp1 = parameters.hbParameters[TsaiAtomicType1][TsaiAtomicType2].valAtp1;
        double valAtp2 = parameters.hbParameters[TsaiAtomicType1][TsaiAtomicType2].valAtp2;

        if (end < 0.1)  // With the current paramters for the solvate energy: end=0 means that hydrogen bonding is not relevent to this pair.
            return null;

        return new HydrogenBondDahiyat(dm, tmpList, p1, p2, end, valAtp1, valAtp2);

    }


}


  
