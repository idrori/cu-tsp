/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.sideChainModelingSolvate;

import meshi.energy.*;
import meshi.util.mathTools.*;
import meshi.molecularElements.atoms.*;
import meshi.geometry.*;
import meshi.parameters.*;

/**
 * This class is a truncated form of the corresponding class in the "meshi.energy.solvation" package. It was designed
 * solely for accelerating the application SCMOD (concurrent sidechain modeling), and should not be used for other
 * purposes.
 * <p/>
 * Please do not use it!
 */

public final class SideChainSolvateEnergy extends CooperativeEnergyTerm {


    /**
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted isOn creating new
     * instances of the arrays.
     */
    private double[] AtomSumSigmC;
    private double[] AtomSumSigmHB;

    /**
     * These fields are for general use in the class
     */
    private int[] lut;     // The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.
    private int atomListSize;
    private SideChainSolvateParametersList parameters; // The instance of the parameter list object.
    private double simpleHBweight; // The weight given to the HB energy term.
    private Spline1D[] splines; // The splines array, that can associate each atom in the protein with the spline that suits its type.
    private SideChainSolvateHBAngle solvateHBAngle; // We create a single instance of SolvateHBAngle, and use it to calculate values, regarding the angular properties of the hydrogen bonds.
    private Atom[] baseAtom; // This array stores pointers to the base getAtoms of every non-hydrogen polar atom in the protein.
    private boolean[] inactive;
    private boolean[] inactiveAtConstructor;
    private double hbEnergy;
    private final Sigma sigma = new Sigma();

    public SideChainSolvateEnergy() {
    }

    /**
     * The constructor parameters:
     * atomList,dm,weight - standard CooperativeEnergyTerm inputs. The 'weight' parameter is the parmeter of
     * the cooperative part (Wcooperative). Note, that changing it will not affect the weight given to the HB
     * part.
     * simpleHBweight - The weight of the HB part.
     * sigmoidBeginsWithH,sigmoidEndsWithH - The transition angles (in degrees) for the HB sigmoid when one of the
     * base getAtoms (or both) is a hydrogen. Above sigmoidEndsWithH the sigmoid is given a value of 1.0 . Bellow
     * sigmoidBeginsWithH the sigmoid is given a value of 0.0 . In between it raises smoothly by cubic spline.
     * sigmoidBeginsNoH,sigmoidEndsNoH - The same as above, only for HB sigmoids where none of the base getAtoms
     * is a hydrogen.
     */
    public SideChainSolvateEnergy(AtomList atomList,
                                  DistanceMatrix dm,
                                  SideChainSolvateParametersList parameters,
                                  double simpleHBweight,
                                  double sigmoidBeginsWithH,
                                  double sigmoidEndsWithH,
                                  double sigmoidBeginsNoH,
                                  double sigmoidEndsNoH,
                                  EnergyInfoElement info) {
        super(toArray(), atomList, dm, parameters, info);
        int c;
        int maxAtomNum = -1;
        comment = "Undefined Solvation";
        atomListSize = atomList.size();
        this.simpleHBweight = simpleHBweight;
        this.parameters = parameters;
        if (parameters == null)
            throw new RuntimeException("The parameters object for this Solvatation term is NULL");
        if (parameters.maxEnd > dm.rMax())
            throw new RuntimeException("This solvatation term can only work if the rMax in " +
                    "the distance matrix is larger than:" + parameters.maxEnd);

        // Creating the auxilary arrays
        AtomSumSigmC = new double[atomListSize];
        AtomSumSigmHB = new double[atomListSize];


        // Creating the lookup table for the atom numbers.
        // The table converts the atom internal number (field of Atom) to its index
        // in the atom list given to the constructor.
        for (c = 0; c < atomListSize; c++) {
            if (atomList.atomAt(c).number() > maxAtomNum)
                maxAtomNum = atomList.atomAt(c).number();
        }
        lut = new int[maxAtomNum + 1];
        for (c = 0; c < maxAtomNum; c++) {
            lut[c] = -1;
        }
        for (c = 0; c < atomListSize; c++) {
            lut[atomList.atomAt(c).number()] = c;
        }

        // Setting up the splines array, that can associate each atom in the protein with the spline that suits its type.
        splines = new Spline1D[atomListSize];
        for (c = 0; c < atomListSize; c++)
            splines[c] = parameters.atomTypeSplines[atomList.atomAt(c).type().ordinal()];

        // Setting up the HB angle object, and the base atom array.
        // See the documentaion of "SolvateHBAngle.java" to see what are the parameters to its
        // constructor.
        solvateHBAngle = new SideChainSolvateHBAngle(dm,
                sigmoidBeginsWithH,
                sigmoidEndsWithH,
                sigmoidBeginsNoH,
                sigmoidEndsNoH);

        // Setting the base atom for every polar non-hydrogen atom. This base atom is used to determined the HB angles.
        // The places in 'baseAtom' that correspond to carbon or hydrogen getAtoms remain 'null'.
        baseAtom = new Atom[atomListSize];
        setBaseAtom();

        // Dealing with inactivation.
        inactive = new boolean[atomListSize];
        inactiveAtConstructor = new boolean[atomListSize];
        for (c = 0; c < atomListSize; c++)
            if (atomList.atomAt(c).name().equals("O") || atomList.atomAt(c).name().equals("N"))
                inactive[c] = false;
            else
                inactive[c] = true;
        evaluate(false, 0, 0);
        for (c = 0; c < atomListSize; c++)
            if ((atomList.atomAt(c).name().equals("O") || atomList.atomAt(c).name().equals("N"))
                    && (AtomSumSigmHB[c] > 0.8))
                inactiveAtConstructor[c] = inactive[c] = true;
            else
                inactiveAtConstructor[c] = inactive[c] = false;
    } // of the constructor

    public void setComment(String str) {
        comment = str;
    }

    public void setSimpleHBweight(double newWeight) {
        simpleHBweight = newWeight;
    }


    public double getNonWeightedHBenergy() {
        return hbEnergy;
    }

    public void inactivateFarFromAtom(Atom atom, double R) {
        double x = atom.x();
        double y = atom.y();
        double z = atom.z();
        Atom atom1;
        for (int c = 0; c < atomListSize; c++) {
            atom1 = atomList.atomAt(c);
            if (!atom1.nowhere()) {
                if (((atom1.x() - x) > R) || ((atom1.x() - x) < -R) ||
                        ((atom1.y() - y) > R) || ((atom1.y() - y) < -R) ||
                        ((atom1.z() - z) > R) || ((atom1.z() - z) < -R))
                    inactive[c] = true;
                else
                    inactive[c] = inactiveAtConstructor[c];
            }
        }
    }


    public void evaluateAtoms() {
        evaluate(true, weight, simpleHBweight);
    }

    public EnergyInfoElement evaluate() {
        return evaluate(false, weight, simpleHBweight);
    }

    public final EnergyInfoElement evaluate(boolean updateAtoms, double cooperativeWeight, double simpleHBweight) {
        double energy = 0;
        double envior;
        double atomEnergy = 0;
        int cc;
        DistanceLists dislist = dm.nonBondedList();
        SideChainSolvateDistanceAttribute sigmaValues;
        Atom atom;
        int ind1, ind2, ind3, ind4;


        //Reseting the auxilary arrays and variables
        hbEnergy = 0.0;
        for (cc = 0; cc < atomListSize; cc++) {
            AtomSumSigmC[cc] = 0;
            AtomSumSigmHB[cc] = 0;
        }

        // First pass over the non-bonded list
        for (DistanceList distanceRow : dislist) {
            for (Distance dis : distanceRow) {
                if (dis.getAttribute(SideChainSolvateDistanceAttribute.SIDE_CHAIN_SOLVATE_ALL_ATOM_ATTRIBUTE) == null) {
                    sigmaValues = new SideChainSolvateDistanceAttribute();
                    dis.addAttribute(sigmaValues);
                } else
                    sigmaValues = (SideChainSolvateDistanceAttribute) dis.getAttribute(SideChainSolvateDistanceAttribute.SIDE_CHAIN_SOLVATE_ALL_ATOM_ATTRIBUTE);
                ind1 = lut[dis.atom1().number()];
                ind2 = lut[dis.atom2().number()];
                updateSigmVals(dis);
                AtomSumSigmC[ind1] += sigmaValues.sigmCa1;
                AtomSumSigmC[ind2] += sigmaValues.sigmCa2;
                AtomSumSigmHB[ind1] += sigmaValues.sigmHBa1;
                AtomSumSigmHB[ind2] += sigmaValues.sigmHBa2;
            }
        }

        //Calculating the energy values. Looping isOn all the getAtoms in the protein.
        for (cc = 0; cc < atomListSize; cc++)
            if (!inactive[cc]) {
                envior = AtomSumSigmC[cc] * (1 - AtomSumSigmHB[cc]); // The functional form of the environment index
                if (envior < 0.0)
                    envior = 0.0;
                splines[cc].calc(envior);
                atomEnergy = cooperativeWeight * splines[cc].s - simpleHBweight * AtomSumSigmHB[cc]; // The energy associated with the atom.
                hbEnergy -= AtomSumSigmHB[cc];
                energy += atomEnergy;
                if (updateAtoms)
                    atomList.atomAt(cc).addEnergy(atomEnergy);
            }

        info.setValue(energy);
        return info;
    }


    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and
     * atom2 in the Distance - dis. The results are updated in the fields of the
     * SolvateDistanceAttribute of dis - sigmaValues.
     */
    private final void updateSigmVals(Distance dis) {
        SideChainSolvateDistanceAttribute sigmaValues =
                (SideChainSolvateDistanceAttribute) dis.getAttribute(SideChainSolvateDistanceAttribute.SIDE_CHAIN_SOLVATE_ALL_ATOM_ATTRIBUTE);
        int TsaiAtomicType1 = parameters.atomicTypeConverter[dis.atom1().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        int TsaiAtomicType2 = parameters.atomicTypeConverter[dis.atom2().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        double tmpSigmaHBatom1, tmpSigma_TagHBatom1, tmpSigmaHBatom2, tmpSigma_TagHBatom2;

        sigmaValues.resetAllSigmVals();

        // Hydrogens are not treated currently
        // -------------------------------
        if (dis.atom1().type().isHydrogen() || dis.atom2().type().isHydrogen()) {
            return;
        }

        // Disactivated getAtoms are not treated
        // -------------------------------
        if (inactive[lut[dis.atom1().number()]] || inactive[lut[dis.atom2().number()]]) {
            return;
        }


        // Calculating the carbon sigmoid of atom1.
        // ---------------------------------------
        if (dis.atom2().type().isCarbon()) {
            sigma.sigma(dis.distance(),
                    parameters.Cend[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.Cp1[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.Cp2[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.CvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.CvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
            sigmaValues.sigmCa1 = sigma.s();
        }

        // Calculating the carbon sigmoid of atom2.
        // ---------------------------------------
        if (dis.atom1().type().isCarbon()) {
            sigma.sigma(dis.distance(),
                    parameters.Cend[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.Cp1[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.Cp2[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.CvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.CvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
            sigmaValues.sigmCa2 = sigma.s();
        }

        // Possible HB Ahoy!!
        // (atom1 and atom2 are capable of forming a hydrogen bond)
        // ------------------
        if ((dis.atom1().type().isOxygen() || dis.atom1().type().isNitrogen() || dis.atom1().type().isSulfur()) &&
                (dis.atom2().type().isOxygen() || dis.atom2().type().isNitrogen() || dis.atom2().type().isSulfur())) {
            // Calculating the HB sigmoid of atom1
            sigma.sigma(dis.distance(),
                    parameters.HBend[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.HBp1[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.HBp2[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.HBvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.HBvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
            tmpSigmaHBatom1 = sigma.s();
            tmpSigma_TagHBatom1 = sigma.s_tag();
            // Calculating the HB sigmoid of atom2
            sigma.sigma(dis.distance(),
                    parameters.HBend[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.HBp1[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.HBp2[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.HBvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.HBvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
            tmpSigmaHBatom2 = sigma.s();
            tmpSigma_TagHBatom2 = sigma.s_tag();
            if ((tmpSigmaHBatom2 > 0) || (tmpSigmaHBatom1 > 0)) {    // Some sort of a HB interaction
                // case 1:   base---O...O---base    or    base---O...N---base
                if (!baseAtom[lut[dis.atom1().number()]].type().isHydrogen() &&
                        !baseAtom[lut[dis.atom2().number()]].type().isHydrogen()) {
                    solvateHBAngle.updateAndEvaluateAtoms1234(
                            baseAtom[lut[dis.atom1().number()]],
                            dis.atom1(),
                            dis.atom2(),
                            baseAtom[lut[dis.atom2().number()]]);
                    sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
                    sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
                } else  // case 2:   N---H...O---base
                    if (baseAtom[lut[dis.atom1().number()]].type().isHydrogen() &&
                            !baseAtom[lut[dis.atom2().number()]].type().isHydrogen()) {
                        solvateHBAngle.updateAndEvaluateAtoms1234(
                                dis.atom1(),
                                baseAtom[lut[dis.atom1().number()]],
                                dis.atom2(),
                                baseAtom[lut[dis.atom2().number()]]);
                        sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
                        sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
                    } else // case 3:   base---O...H---N
                        if (!baseAtom[lut[dis.atom1().number()]].type().isHydrogen() &&
                                baseAtom[lut[dis.atom2().number()]].type().isHydrogen()) {
                            solvateHBAngle.updateAndEvaluateAtoms1234(
                                    baseAtom[lut[dis.atom1().number()]],
                                    dis.atom1(),
                                    baseAtom[lut[dis.atom2().number()]],
                                    dis.atom2());
                            sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
                            sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
                        } else // case 4:   N---H...H---N
                            if (baseAtom[lut[dis.atom1().number()]].type().isHydrogen() &&
                                    baseAtom[lut[dis.atom2().number()]].type().isHydrogen()) {
                                solvateHBAngle.updateAndEvaluateAtoms1234(
                                        dis.atom1(),
                                        baseAtom[lut[dis.atom1().number()]],
                                        baseAtom[lut[dis.atom2().number()]],
                                        dis.atom2());
                                sigmaValues.sigmHBa1 = tmpSigmaHBatom1 * solvateHBAngle.hbAngScore();
                                sigmaValues.sigmHBa2 = tmpSigmaHBatom2 * solvateHBAngle.hbAngScore();
                            }
            } // Of checkin if a HB interaction soes exist
        } // Of checking for HB forming pair
    } // of updateSigmVals


    /**
     * For any polar atom (O,N or S in Cys) we calculate the base atom that participate in the definition of the hydrogen
     * bond angle. This base atom is the attached hydrogen (if present), or the heavy atom to which the polar atom is
     * attached (when the hydrogen is not present).
     */
    private final void setBaseAtom() {
        Atom atom1, atom2;
        for (int c1 = 0; c1 < atomList.size(); c1++) {
            atom1 = atomList.atomAt(c1);
            for (int c2 = 0; c2 < atom1.bonded().size(); c2++) {
                atom2 = atom1.bonded().atomAt(c2);
                // Treating OXYGEN getAtoms.
                // This is easy because the oxygens in proteins
                // are always tied to one atom only.
                if (atom1.type().isOxygen())
                    baseAtom[lut[atom1.number()]] = atom2;
                // Treating NITROGEN getAtoms.
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
                // Treating the SULFUR getAtoms of Cystines.
                if (atom1.type() == AtomType.CSG)
                    baseAtom[lut[atom1.number()]] = atom2;
            }
        }
    }

}
	
