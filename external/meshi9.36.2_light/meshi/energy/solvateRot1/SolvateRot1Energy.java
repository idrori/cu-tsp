/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvateRot1;

import meshi.energy.CooperativeEnergyTerm;
import meshi.energy.EnergyInfoElement;
import meshi.energy.solvation.SigmaParameters;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.mathTools.Spline1D;
import meshi.util.rotamericTools.Rot1Arrays;
import meshi.util.rotamericTools.RotamericTools;

import java.util.Iterator;


public final class SolvateRot1Energy extends CooperativeEnergyTerm implements Rot1Arrays {


    /**
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted isOn creating new
     * instances of the arrays.
     */
    private double[] residueVals;
    private double[] splineVals;
    private double[] dSplineVals;
    private double[] AtomSumSigmC;
    private double[] dAtomEnergydEnvior;
    private double[] forceX;
    private double[] forceY;
    private double[] forceZ;
    /**
     * These fields are for general use in the class
     */
    private int[] lut;     // The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.
    private int atomListSize;
    private SolvateRot1ParametersList parameters; // The instance of the parameter list object.
    private Spline1D[] splines; // The splines array, that can associate each atom in the protein with the spline that suits its type.
    private double[][] means;
    private double[][] oneOverSTDs;
    private double maximalZ;


    public SolvateRot1Energy(AtomList atomList,
                             DistanceMatrix dm,
                             SolvateRot1ParametersList parameters,
                             double[][] pp,
                             double maximalZ,
                             EnergyInfoElement info) {
        super(toArray(dm), atomList, dm, parameters, info);
        int c;
        int maxAtomNum = -1;
        comment = "Undefined Solvation";
        atomListSize = atomList.size();
        this.parameters = parameters;
        this.maximalZ = maximalZ;
        if (parameters == null)
            throw new RuntimeException("The parameters object for this Solvatation term is NULL");
        if (parameters.maxEnd > dm.rMax())
            throw new RuntimeException("This solvatation term can only work if the rMax in " +
                    "the distance matrix is larger than:" + parameters.maxEnd);
        calcMeansAndStd(pp);
        residueVals = new double[pp.length];

        // Creating the auxilary arrays
        splineVals = new double[atomListSize];
        dSplineVals = new double[atomListSize];
        AtomSumSigmC = new double[atomListSize];
        dAtomEnergydEnvior = new double[atomListSize];
        forceX = new double[atomListSize];
        forceY = new double[atomListSize];
        forceZ = new double[atomListSize];


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

    } // of the constructor

    public void setComment(String str) {
        comment = str;
    }

    public void evaluateAtoms() {
        evaluate(true);
    }

    public EnergyInfoElement evaluate() {
        return evaluate(false);
    }

    public final EnergyInfoElement evaluate(boolean updateAtoms) {
        double energy = 0;
        double atomEnergy = 0;
        int iAtom;
        DistanceLists nonBondedList = dm.nonBondedList();
        Iterator iter;
        SolvationRot1DistanceAttribute sigmaValues;
        Atom atom;
        int ind1, ind2;


        //Reseting the auxilary arrays and variables
        for (iAtom = 0; iAtom < atomListSize; iAtom++) {
            AtomSumSigmC[iAtom] = 0;
            forceX[iAtom] = forceY[iAtom] = forceZ[iAtom] = 0.0;
        }
        for (iAtom = 0; iAtom < residueVals.length; iAtom++) {
            residueVals[iAtom] = 0.0;
        }

        // First pass over the non-bonded list
        iter = nonBondedList.iterator();
        for (DistanceList row : nonBondedList) {
            for (Distance distance : row) {
                if (distance.getAttribute(SolvationRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE) == null) {
                    sigmaValues = createSigmaValues(distance);
                    distance.addAttribute(sigmaValues);
                    if (distance.mode().frozen)
                        updateSigmVals(distance);
                }
                if (!distance.mode().frozen) {
                    updateSigmVals(distance);
                    sigmaValues = (SolvationRot1DistanceAttribute) distance.getAttribute(SolvationRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);
                    ind1 = lut[distance.atom1().number()];
                    ind2 = lut[distance.atom2().number()];
                    AtomSumSigmC[ind1] += sigmaValues.sigmCa1;
                    AtomSumSigmC[ind2] += sigmaValues.sigmCa2;
                }
            }
        }

        //Calculating the energy values. Looping isOn all the getAtoms in the protein... TWICE!
        for (iAtom = 0; iAtom < atomListSize; iAtom++) {
            if (!atomList.atomAt(iAtom).frozen()) {
                splines[iAtom].calc(AtomSumSigmC[iAtom]);
                splineVals[iAtom] = splines[iAtom].s;
                dSplineVals[iAtom] = splines[iAtom].s_tag;
                residueVals[atomList.atomAt(iAtom).residueNumber()] += splines[iAtom].s;
            }
        }
        // A second time
        for (iAtom = 0; iAtom < atomListSize; iAtom++) {
            atom = atomList.atomAt(iAtom);
            if (!atom.frozen()) { // This line was deleted isOn 19.6.2006. It looks like a bug...
                if (means[atom.residueNumber()] != null) {
                    if (atom.name.equals("CA")) {
                        energy += weight *
                                Math.min((residueVals[atom.residueNumber()] -
                                        means[atom.residueNumber()][0]) *
                                        oneOverSTDs[atom.residueNumber()][0], maximalZ);
                        if (updateAtoms)
                            atom.addEnergy(weight *
                                    Math.min((residueVals[atom.residueNumber()] -
                                            means[atom.residueNumber()][0]) *
                                            oneOverSTDs[atom.residueNumber()][0], maximalZ));
                    }
                    dAtomEnergydEnvior[iAtom] = weight * oneOverSTDs[atom.residueNumber()][0] * dSplineVals[iAtom];
                }
            } // This line was deleted isOn 19.6.2006. It looks like a bug...
        }

        // Second pass over the non-bonded list
        for (DistanceList row : nonBondedList) {
            for (Distance distance : row) {
                if (!distance.mode().frozen) {
                    sigmaValues = (SolvationRot1DistanceAttribute) distance.getAttribute(SolvationRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);
                    ind1 = lut[distance.atom1().number()];
                    ind2 = lut[distance.atom2().number()];

                    // Doing the self derivatives
                    forceX[ind1] += dAtomEnergydEnvior[ind1] * sigmaValues.dsigmCa1dx1;
                    forceY[ind1] += dAtomEnergydEnvior[ind1] * sigmaValues.dsigmCa1dy1;
                    forceZ[ind1] += dAtomEnergydEnvior[ind1] * sigmaValues.dsigmCa1dz1;

                    forceX[ind2] += dAtomEnergydEnvior[ind2] * sigmaValues.dsigmCa2dx2;
                    forceY[ind2] += dAtomEnergydEnvior[ind2] * sigmaValues.dsigmCa2dy2;
                    forceZ[ind2] += dAtomEnergydEnvior[ind2] * sigmaValues.dsigmCa2dz2;

                    // Doing the cross derivatives
                    forceX[ind2] += dAtomEnergydEnvior[ind1] * sigmaValues.dsigmCa1dx2;
                    forceY[ind2] += dAtomEnergydEnvior[ind1] * sigmaValues.dsigmCa1dy2;
                    forceZ[ind2] += dAtomEnergydEnvior[ind1] * sigmaValues.dsigmCa1dz2;

                    forceX[ind1] += dAtomEnergydEnvior[ind2] * sigmaValues.dsigmCa2dx1;
                    forceY[ind1] += dAtomEnergydEnvior[ind2] * sigmaValues.dsigmCa2dy1;
                    forceZ[ind1] += dAtomEnergydEnvior[ind2] * sigmaValues.dsigmCa2dz1;
                }
            }
        }// second pass isOn the non bonded list

        // Finally, the appropriate forces are assigned for every atom.
        for (iAtom = 0; iAtom < atomListSize; iAtom++) {
            atom = atomList.atomAt(iAtom);
            if (!atom.frozen()) {
                atom.addToFx(-forceX[iAtom]); // Negating so that it is realy force (and not a mere derivative)
                atom.addToFy(-forceY[iAtom]); // Negating so that it is realy force
                atom.addToFz(-forceZ[iAtom]); // Negating so that it is realy force
            }
        }
        info.setValue(energy);
        return info;
    }


    /**
     * Here goes the processing of the various sigmoid functions concerning atom1 and
     * atom2 in the Distance - dis. The results are updated in the fields of the
     * SolvateDistanceAttribute of dis - sigmaValues.
     */
    private SolvationRot1DistanceAttribute createSigmaValues(Distance distance) {
        if (distance.atom1().type().isHydrogen() || distance.atom2().type().isHydrogen()) return null;

        int TsaiAtomicType1 = parameters.atomicTypeConverter[distance.atom1().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        int TsaiAtomicType2 = parameters.atomicTypeConverter[distance.atom2().type().ordinal()]; // Converting from the 190 atom types to the 14 defined in Tsai 99'
        SigmaParameters sigmaParameters_12, sigmaParameters_21;

        if (distance.atom2().type().isCarbon()) {
            sigmaParameters_12 = new SigmaParameters(parameters.Cend[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.Cp1[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.Cp2[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.CvalAtp1[TsaiAtomicType1][TsaiAtomicType2],
                    parameters.CvalAtp2[TsaiAtomicType1][TsaiAtomicType2]);
        } else sigmaParameters_12 = null;
        if (distance.atom1().type().isCarbon()) {
            sigmaParameters_21 = new SigmaParameters(parameters.Cend[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.Cp1[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.Cp2[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.CvalAtp1[TsaiAtomicType2][TsaiAtomicType1],
                    parameters.CvalAtp2[TsaiAtomicType2][TsaiAtomicType1]);
        } else sigmaParameters_21 = null;
        return new SolvationRot1DistanceAttribute(distance,sigmaParameters_12,sigmaParameters_21);
    }
        private final void updateSigmVals(Distance dis) {
        SolvationRot1DistanceAttribute sigmaValues =
                (SolvationRot1DistanceAttribute) dis.getAttribute(SolvationRot1DistanceAttribute.SOLVATE_ROT1_ATTRIBUTE);

        sigmaValues.resetAllSigmVals();

        // Hydrogens are not treated currently
        // -------------------------------
        if (dis.atom1().type().isHydrogen() || dis.atom2().type().isHydrogen()) {
            return;
        }


        // Calculating the carbon sigmoid of atom1.
        // ---------------------------------------
        if (dis.atom2().type().isCarbon()) {
            sigmaValues.sigma_12.sigma(dis.distance());
            sigmaValues.sigmCa1 = sigmaValues.sigma_12.s();
            sigmaValues.dsigmCa1dx1 = sigmaValues.sigma_12.s_tag() * dis.dDistanceDx();
            sigmaValues.dsigmCa1dy1 = sigmaValues.sigma_12.s_tag() * dis.dDistanceDy();
            sigmaValues.dsigmCa1dz1 = sigmaValues.sigma_12.s_tag() * dis.dDistanceDz();
            sigmaValues.dsigmCa1dx2 = -sigmaValues.dsigmCa1dx1;
            sigmaValues.dsigmCa1dy2 = -sigmaValues.dsigmCa1dy1;
            sigmaValues.dsigmCa1dz2 = -sigmaValues.dsigmCa1dz1;
        }

        // Calculating the carbon sigmoid of atom2.
        // ---------------------------------------
        if (dis.atom1().type().isCarbon()) {
            sigmaValues.sigma_21.sigma(dis.distance());
            sigmaValues.sigmCa2 = sigmaValues.sigma_21.s();
            sigmaValues.dsigmCa2dx1 = sigmaValues.sigma_21.s_tag() * dis.dDistanceDx();
            sigmaValues.dsigmCa2dy1 = sigmaValues.sigma_21.s_tag() * dis.dDistanceDy();
            sigmaValues.dsigmCa2dz1 = sigmaValues.sigma_21.s_tag() * dis.dDistanceDz();
            sigmaValues.dsigmCa2dx2 = -sigmaValues.dsigmCa2dx1;
            sigmaValues.dsigmCa2dy2 = -sigmaValues.dsigmCa2dy1;
            sigmaValues.dsigmCa2dz2 = -sigmaValues.dsigmCa2dz1;
        }

    } // of updateSigmVals


    public void calcMeansAndStd(double[][] pp) {
        means = new double[pp.length][];
        oneOverSTDs = new double[pp.length][];
        for (int c = 0; c < pp.length; c++) {
            if ((pp[c] != null) && (((int) pp[c][2]) != 0) && (((int) pp[c][2]) != 5)) {
                double[] tmp = RotamericTools.getMean(((int) pp[c][2]), pp[c][3], 0);
                means[c] = new double[1];
                oneOverSTDs[c] = new double[1];
                means[c][0] = tmp[0];
                oneOverSTDs[c][0] = 1 / (tmp[1] + 0.00000000001);
            }
        }
    }

}
	
