/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation;

import meshi.energy.CooperativeEnergyTerm;
import meshi.energy.EnergyInfoElement;
import meshi.energy.TotalEnergy;
import meshi.energy.solvation.hydrogenBonds.AbstractHydrogenBondList;
import meshi.energy.solvation.hydrogenBonds.HydrogenBondDahiyatList;
import meshi.geometry.*;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;
import meshi.parameters.ResidueType;
import meshi.util.Stat;
import meshi.util.Utils;
import meshi.util.mathTools.Spline2D;

/**
 * The implementation of the cooperative solvation term for proteins as described in Kalisman & Keasar (2008).
 * Since the cooperative solvation is described in the above paper, we bring here only the implementations
 * details. Especially regarding the calculation of the derivatives, which was too lengthy for the paper.
 * <p/>
 * The class allows to give a different weight to the solvation energies of certain atom types in the final
 * summation described in Eq. 5. These atom types are side-chain polars, side-chain carbons, and backbone
 * polars. The class also include a regular hydrogen bond term. The functional form is therefore:
 * <p/>
 * Esolv = weightSCPolarSolvation*Eside_chain_polars +
 * weightSCCarbonSolvation*Eside_chain_carbons +
 * weightBBPolarSolvation*Ebackbone_polars +
 * weightHB*Ehb
 * <p/>
 * Where Ehb is the negative of the HBC summation over all getAtoms. The weights are defined in the "Creator" class,
 * and passed to the constructor as parameters.
 * <p/>
 * General remarks:
 * ----------------
 * 1) This term is derivable twice (because the splines are derivable twice).
 * 2) This term is using hydrogen bond description that is dependent isOn two angles in the bond. This decription
 * follows that of McDonald and Thornton (1994) and Dahiyat et al. (1996). The only place where the hydrogen bond
 * list is declared explicitly is in line 204. This means that any hydrogen bond implementation that extends the
 * "AbstractHydrogenBondList" template can be used, by correcting line 204.
 * 3) We calculate the regular hydrogen bond energy term (Ehb) together with the solvation terms themselves, since
 * the hydrogen bonds calculation is a by-product of the first step in the solvation evaluation, and is thus for free.
 * 4) Disulfide bonds are treated as "hydrogen bonds" between the SG's of two cystines.
 * 5) The SD sulfor of methionine is treated as a hydrophobic carbon.
 * 6) See the remarks in the "Creators" classes for a quick start isOn how to create a working instance of the term.
 * <p/>
 * <p/>
 * The energy evaluation:
 * ----------------------
 * The energy value and derivatives is calculated in 3 steps:
 * 1) A first pass over the non-bonded list. Each Distance instance in the non-bonded-list, is used to update the
 * CNC's and HBC's of its atom pairs (Eqs. 1 and 2, respectively). The partial derivatives of the CNC's and HBC's
 * with respect to the distance getAtoms are alsoclaculated.  Since some of this values will be needed also in step 3,
 * we save them in an instance of "SolvateDistanceAttribute" that is attached as an "Attribute" to the Distance instance.
 * 2) A pass isOn the atom list. Once we have the CNC and HBC of every atom in the protein, we can proceed to calculate the solvation energy
 * associated with every atom. In this implementation we combined the EI(CNC,HBC) evaluation (Eq. 3) of every atom and
 * the -log(spline(EI)) evaluation (Eq. 4) into a single step by using a 2D spline, i.e. spline2D(CNC,HBC). The 2D spline
 * is, of course, atom type specific. The derivatives of each atom solvation energy value with respect to the HBC and CNC
 * are also calculated.
 * 3) A second pass over the non-bonded list. Equiped with the The derivatives of the atom energies (with respect
 * to the CNC's and HBC's) from step 2, we can now calculate the energy derivative with respect to the atomic
 * coordinates. In this step we simply make sure that every term that arises from the derivative chain rule is accounted for.
 */
public class SolvationEnergy extends CooperativeEnergyTerm {

    public static final double BETA = 1;
    public static boolean debug = false;
    // Relative strength of SALT BRIDGES compared with regular HYDROGEN BONDS for desolvation purposes, i.e. in
    // regard to the effect isOn observed CNC medians. Following Table 1 in the paper.
    public static final double SALT_BRIDGE_STRENGTH_ASP_OD = 0.45;
    public static final double SALT_BRIDGE_STRENGTH_GLU_OE = 0.45;
    public static final double SALT_BRIDGE_STRENGTH_LYS_NZ = 1.65;
    public static final double SALT_BRIDGE_STRENGTH_ARG_NH = 1.5;
    public static final double SALT_BRIDGE_STRENGTH_TRO = 0.45;
    public static final double SALT_BRIDGE_STRENGTH_TRN = 1.8;
    // A flag for distances that do not take part in this term
    public static final SolvationDistanceAttribute NOT_PATRICIPATING_DISTANCE = new SolvationDistanceAttribute(null,null);
    /**
     * The following parameter allow for a different weighting of SALT BRIDGES compared with regular HYDROGEN BONDS for the
     * Ehb energy, that is also claculated.
     */
    public static final double SALT_BRIDGE_STRENGTH_GENERAL = 1.0;

    /**
     * These are fields for temporary array results that are needed in the evaluation stage.
     * They are declared as fields so that time will not be waisted isOn creating new
     * instances of the arrays. The lengths of these arrays is the length of the atom list.
     * The indexing to these arrays is by the index of the atom in the atom list.
     */
    private double[] CNC;  // Eq. 1
    public double[] cnc() {return CNC;}
    private double[] HBC;  // Eq. 2
    public double[] hbc() {return HBC;}
    public double[][] dSpline;
    public double[] dStdDe;
    public double[][] forceXYZ;
    public PerAtomWeights[] perAtomWeights;

    private double energy = 0;
    private double energySCPolarTotal = 0;
    private double energySCCarbonTotal = 0;
    private double energyBBPolarTotal = 0;
    private double energyBBCarbonTotal = 0;
    private double energyHbTotal = 0;
    private double energyBuriedHbTotal = 0;
    //For statistics
    private int nBbPolar, nBbCarbon, nScPolar, nScCarbon;

    protected final SolvationDistanceAttributeCreator solvationDistanceAttributeCreator   = new SolvationDistanceAttributeCreator();
    protected final SolvationDistanceAttributeCreator solvationDistanceAttributeCreatorSS = new SolvationDistanceAttributeCreatorSS();
    public SolvationDistanceAttributeCreator solvationDistanceAttributeCreator() {return solvationDistanceAttributeCreator;}
    /**
     * These following fields are for general use in the class
     **/
    /** Size of the atom list. **/
    /**
     * The instance of the parameter list object. *
     */
    private SolvateParametersList parameters;
    /**
     * The look-up table (lut) converts the atom internal number (field of Atom), which is the index of the array, to its
     * index in the atom list given to the constructor. *
     */
    private int[] lut;

    /**
     * Setting the general type for each atom in the atom list: (0) Carbon (1) Backbone polar, (2) Sidechain polar (3) Hydrogens
     * This is done to save time isOn type checking. The index to the array is the index of an atom in the atom list. *
     */
    private enum SuperType {
        BB_CARBON, //A backbone carbon atom.
        BB_POLAR,
        SC_POLAR,
        SC_CARBON,
        HYDROGEN,
        OTHER;


        public static SuperType type(Atom atom) {
            AtomType atomType = atom.type();
            ResidueType residueType = atom.residue().type();
            String atomName = atom.name();
            if (atom.isHydrogen()) return HYDROGEN;
            if ((atom.type() == AtomType.TRN) || (atom.type() == AtomType.TRO) || (atom.type() == AtomType.TRC))
                return OTHER;

            if (atomType.isCarbon() || (atomType == AtomType.MSD)) {
                if (atomType.backbone()) return BB_CARBON;
                else return SC_CARBON;
            } else {
                if (atomType.backbone()) return BB_POLAR;
                else return SC_POLAR;
            }
        }
    }

    /**
     * The 2D spline array (i.e. spline(CNC,HBC)) The indexing in the array is the
     * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) *
     */
    private Spline2D[] splines;
    public Spline2D[] splines() {return splines;}

    /**
     * The only hydrogen bond list class in the term. *
     */
    public AbstractHydrogenBondList solvateHB;


    // Weights
    private double weightSCPolarSolvation;
    private double weightBBPolarSolvation;
    private double weightSCCarbonSolvation;
    private double weightBBCarbonSolvation;
    private double weightHbSolvation;
    private double weightStdSolvation;

    public final MolecularSystem molecularSystem;
    public final int molecularSystemSize;
    private Atom atom;
    private int atomTypeNumber;
    private AtomType atomType;
    private Stat stat;

    public SolvationEnergy() {
        molecularSystem = null;
        molecularSystemSize = -1;
    }

    public double[] atomEnergies;



    /**
     * See the comment at the top of the class for descriptions isOn the weights.
     */

    public SolvationEnergy(AtomList atomList,
                           DistanceMatrix dm,
                           SolvateParametersList parameters,
                           EnergyInfoElement info){
        super(toArray(dm), atomList, dm, parameters, info);
         SolvationInfoElement solvationInfo = (SolvationInfoElement) info;

        molecularSystem = atomList.molecularSystem();

        weightSCPolarSolvation     = solvationInfo.scPolarInfo.weight();
        weightBBPolarSolvation     = solvationInfo.bbPolarInfo.weight();
        weightSCCarbonSolvation    = solvationInfo.scCarbonInfo.weight();
        weightBBCarbonSolvation    = solvationInfo.bbCarbonInfo.weight();
        weightHbSolvation          = solvationInfo.hbInfo.weight();
        weightStdSolvation         = solvationInfo.stdInfo.weight();

        this.comment = "Solvation";
        int iAtom;

        molecularSystemSize = molecularSystem.size();
        this.parameters = parameters;
        if (parameters == null)
            throw new RuntimeException("The parameters object for this Solvatation term is NULL");
        if (parameters.maxEnd > dm.rMax())
            throw new RuntimeException("This solvatation term can only work if the rMax in " +
                    "the distance matrix is larger than:" + parameters.maxEnd);

        atomEnergies = new double[molecularSystemSize];

        // Creating the auxilary arrays
        CNC = new double[molecularSystemSize];
        HBC = new double[molecularSystemSize];
        perAtomWeights = new PerAtomWeights[molecularSystemSize];
        dSpline = new double[molecularSystemSize][2];
        forceXYZ = new double[molecularSystemSize][3];
        dStdDe   = new double[molecularSystemSize];

        // Setting up the splines array, that can associate each residue in the protein with the spline that
        // suits its type.
        splines = new Spline2D[molecularSystemSize];
        for (iAtom = 0; iAtom < molecularSystemSize; iAtom++) {
            atom = molecularSystem.get(iAtom).atom;
            if (!atom.nowhere()) {
                atomTypeNumber = atom.type().ordinal();
                splines[iAtom] = parameters.splines[atomTypeNumber];
            }
        }

        // Setting up the HB list
        solvateHB = new HydrogenBondDahiyatList(dm, atomList, parameters);
        updateableResources.add(solvateHB);

        // Determining the superType for each atom: (0) Hydrophobic, (1) Backbone polar, (2) Sidechain polar (3) Hydrogens
        //For statistics
        nBbPolar = 0;
        nBbCarbon = 0;
        nScPolar = 0;
        nScCarbon = 0;
        for (iAtom = 0; iAtom < molecularSystemSize; iAtom++) {
            atom = molecularSystem.get(iAtom).atom;
            SuperType superType = SuperType.type(atom);
            switch (superType) {
                case SC_CARBON:
                    perAtomWeights[iAtom] = new PerAtomWeights(weightSCCarbonSolvation, weightSCCarbonSolvation, 0, 0, 0);
                    nScCarbon++;
                    break;
                case SC_POLAR:
                    perAtomWeights[iAtom] = new PerAtomWeights(weightSCPolarSolvation, 0, weightSCPolarSolvation, 0, 0);
                    nScPolar++;
                    break;
                case BB_POLAR:
                    perAtomWeights[iAtom] = new PerAtomWeights(weightBBPolarSolvation, 0, 0, weightBBPolarSolvation, 0);
                    nBbPolar++;
                    break;
                case BB_CARBON:
                    perAtomWeights[iAtom] = new PerAtomWeights(weightBBCarbonSolvation, 0, 0, 0, weightBBCarbonSolvation);
                    nBbCarbon++;
                    break;
                default:
                    perAtomWeights[iAtom] = null;
            }
            stat = new Stat();
        }
    }  // Of constructor

    public void evaluateAtoms() {
        evaluate(true);
    }

    /**
     * Calculates Esolv with the weights given in the constructor.
     */
    public EnergyInfoElement evaluate() {
        return evaluate(false);
    }

    /**
     * Calculates Esolv with the weights you give as parameters!
     */
    public final EnergyInfoElement evaluate(boolean updateAtoms) {

        double energy;
        stat.reset();
        resetAuxilaryArrays();
        firstPassOverTheNonBondedList();
        energy = calculatingTheSolvationEnergiesOfEachAtom(updateAtoms);
        secondPassOverTheNonBondedList();
        assignForcesToAtoms();
        info.setValue(energy);
        return info;
    }

    private void resetAuxilaryArrays() {
        int iAtom;
        //Reseting the auxilary arrays
        for (iAtom = 0; iAtom < molecularSystem.size(); iAtom++) {
            CNC[iAtom] = 0;
            HBC[iAtom] = 0;
            dSpline[iAtom][0] = 0;
            dSpline[iAtom][1] = 0;
            forceXYZ[iAtom][0] = forceXYZ[iAtom][1] = forceXYZ[iAtom][2] = 0.0;
        }
    }

    public void firstPassOverTheNonBondedList() {
        firstPassOverTheNonBondedList(false);
    }

    public void firstPassOverTheNonBondedList(boolean debug) {

        int atomNumber1, atomNumber2;
        SolvationDistanceType type;
        energyHbTotal = 0;
        energyBuriedHbTotal = 0;

        DistanceLists nonBondedList = dm.nonBondedList();
        SolvationDistanceAttribute sigmaValues;
        SolvationDistanceAttributePolars sigmaValuesP;

        // ***********************************
        // First pass over the non-bonded list
        // ***********************************
        for (DistanceList distanceRow : nonBondedList) {
            if (!distanceRow.atomOneHeavy) continue;
            atomNumber1 = distanceRow.atomOneNumber;
            for (Distance distance : distanceRow) {
                if (!distance.heavy) continue;
                sigmaValues = (SolvationDistanceAttribute) distance.getAttribute(SolvationDistanceAttribute.SOLVATION_ALL_ATOM_ATTRIBUTE);
                if (sigmaValues == null) {
                    sigmaValues = solvationDistanceAttributeCreator().create(distance, parameters, solvateHB);
                    distance.addAttribute(sigmaValues);
                }
                if (sigmaValues == NOT_PATRICIPATING_DISTANCE) continue;
                type = sigmaValues.type;
                sigmaValues.update();
                //Distances with hydrogens are ignored.
                atomNumber2 = distance.atom2.number;

                switch (type) {
                    case HYDROPHOBIC_POLAR:
                        HBC[atomNumber1] += sigmaValues.sigmCa12;
                        CNC[atomNumber2] += sigmaValues.sigmCa21;
                        break;
                    case POLAR_HYDROPHOBIC:
                        CNC[atomNumber1] += sigmaValues.sigmCa12;
                        HBC[atomNumber2] += sigmaValues.sigmCa21;
                        break;
                    case TWO_HYDROPHOBIC:
                        CNC[atomNumber1] += sigmaValues.sigmCa12;
                        CNC[atomNumber2] += sigmaValues.sigmCa21;
                        break;
                    case TWO_POLARS:
                        sigmaValuesP = (SolvationDistanceAttributePolars) sigmaValues;
                        sigmaValuesP.update();
                        energyHbTotal += sigmaValuesP.sigmHB;
                        if (CNC[atomNumber1]>6) energyBuriedHbTotal += 0.5*sigmaValuesP.sigmHB;
                        if (CNC[atomNumber2]>6) energyBuriedHbTotal += 0.5*sigmaValuesP.sigmHB;

                        HBC[atomNumber1] += sigmaValuesP.saltBridgeFactorA1 * sigmaValuesP.sigmHB;
                        HBC[atomNumber2] += sigmaValuesP.saltBridgeFactorA2 * sigmaValuesP.sigmHB;
                        break;
                }
            }
        }
    }

    public void calculateAtomEnergy(AtomCore atom, boolean updateAtoms) {
        int atomNumber;
        if (atom.atom.isHydrogen()) return;
        if (atom.atom.nowhere())    return;
        atomNumber = atom.number;
        Spline2D spline = splines[atomNumber];
        double s;
        double s_tag_x;
        double s_tag_y;
        PerAtomWeights weights;



        if (spline.calc(CNC[atomNumber], HBC[atomNumber])) {
            s = spline.s;
            s_tag_x = spline.s_tag_x;
            s_tag_y = spline.s_tag_y;
        } else {
            s = 99999;
            s_tag_x = 0;
            s_tag_y = 0;
        }
        weights = perAtomWeights[atomNumber];
        if (weights != null) {
            double myEnergy = s * weights.weight;
            energySCPolarTotal += s * weights.sideChainPolarWeight;
            energySCCarbonTotal += s * weights.sideChainCarbonWeight;
            energyBBPolarTotal += s * weights.polarBackboneWeight;
            energyBBCarbonTotal += s * weights.polarSideChainCarbonWeight;
            dSpline[atomNumber][0] = s_tag_x * weights.weight;
            dSpline[atomNumber][1] = s_tag_y * weights.weight;
            stat.add(myEnergy);
            energy += myEnergy; // We want to count each HB once.
            if (updateAtoms)
                molecularSystem.get(atomNumber).atom.addEnergy(myEnergy);
            if (molecularSystem.get(atomNumber).atom.number() != atomNumber)
                throw new RuntimeException("Check the size of atomEnergies");
            atomEnergies[atomNumber] += myEnergy;
        }
    }
    public double calculatingTheSolvationEnergiesOfEachAtom(boolean updateAtoms) {
        energy = 0;
        energySCPolarTotal = 0;
        energySCCarbonTotal = 0;
        energyBBPolarTotal = 0;
        energyBBCarbonTotal = 0;
        double e, sumE = 0, sumE2 = 0, stdE, entropy = 0, sumExpE = 0 ;
        int n = molecularSystemSize;
        double invN = 1.0/n;
        double[] expE = new double[molecularSystem.size()];

        // ***************************************************************************
        // A pass over the atom list. Calculating the solvation energies of each atom.
        // ***************************************************************************
        for (int i = 0; i < atomEnergies.length; i++) atomEnergies[i] = 0;

        for (AtomCore atom : molecularSystem) {
            calculateAtomEnergy(atom, updateAtoms);
            e = atomEnergies[atom.number];
            sumE += e;
            sumE2 += e * e;
            expE[atom.number] = Math.exp(e*BETA);
            sumExpE += expE[atom.number];
        }
        for (AtomCore atom : molecularSystem) {
            if (!atom.atom.isHydrogen()) {
                double p = expE[atom.number] / sumExpE;
                entropy -= p * Math.log(p);
            }
        }
        stdE = Math.sqrt(sumE2*invN-(sumE*invN)*(sumE*invN));
        energy += stdE*weightStdSolvation;
        for (AtomCore atom : molecularSystem) {
            e = atomEnergies[atom.number];
            dStdDe[atom.number] = (0.5 / stdE) * 2 * invN * (e - invN * sumE)*weightStdSolvation;
        }
        ((SolvationInfoElement)info).scPolarInfo.setValue(energySCPolarTotal);
        ((SolvationInfoElement)info).scCarbonInfo.setValue(energySCCarbonTotal);
        ((SolvationInfoElement)info).bbPolarInfo.setValue(energyBBPolarTotal);
        ((SolvationInfoElement)info).bbCarbonInfo.setValue(energyBBCarbonTotal);
        ((SolvationInfoElement)info).hbInfo.setValue(energyHbTotal);
        ((SolvationInfoElement)info).buriedHbinfo.setValue(energyBuriedHbTotal);
        ((SolvationInfoElement)info).stdInfo.setValue(stat.getStd());
        ((SolvationInfoElement)info).entropyInfo.setValue(entropy);
//For statistics
        ((SolvationInfoElement)info).bbPolarN.setValue(new Double(nBbPolar));
        ((SolvationInfoElement)info).bbCarbonN.setValue(new Double(nBbCarbon));
        ((SolvationInfoElement)info).scPolarN.setValue(new Double(nScPolar));
        ((SolvationInfoElement)info).scCarbonN.setValue(new Double(nScCarbon));
         return energy;
    }

    public void secondPassOverTheNonBondedList() {
        int atomNumber1, atomNumber2;
        SolvationDistanceType type;
        DistanceLists nonBondedList = dm.nonBondedList();
        SolvationDistanceAttribute sigmaValues = null;
        SolvationDistanceAttributePolars sigmaValuesBP = null;
        // ************************************
        // Second pass over the non-bonded list
        // ************************************
        double forceX,forceY,forceZ;
        for (DistanceList distanceRow : nonBondedList) {
            if (!distanceRow.atomOneHeavy) continue;
            atomNumber1 = distanceRow.atomOneNumber;
            double dSplineDCNC1  = dSpline[atomNumber1][0];
            double dSplineDHBC1  = dSpline[atomNumber1][1];
            double dStdDe1       = dStdDe[atomNumber1];
            double[] forceXYZ1   = forceXYZ[atomNumber1];
            for (Distance distance : distanceRow) {
                if (!distance.heavy) continue;
                sigmaValues = (SolvationDistanceAttribute) distance.getAttribute(SolvationDistanceAttribute.SOLVATION_ALL_ATOM_ATTRIBUTE);
                if (sigmaValues == NOT_PATRICIPATING_DISTANCE) continue;
                type = sigmaValues.type;
                atomNumber2 = distance.atom2.number;
                double[] forceXYZ2   = forceXYZ[atomNumber2];
                double dSplineDCNC2  = dSpline[atomNumber2][0];
                double dSplineDHBC2  = dSpline[atomNumber2][1];
                double dStdDe2       = dStdDe[atomNumber2];
                switch (type) {
                    case HYDROPHOBIC_POLAR:
                        forceX = dSplineDHBC1 * sigmaValues.dsigmCa12dx * (1 + dStdDe1);
                        forceY = dSplineDHBC1 * sigmaValues.dsigmCa12dy * (1 + dStdDe1);
                        forceZ = dSplineDHBC1 * sigmaValues.dsigmCa12dz * (1 + dStdDe1);
                        forceX += dSplineDCNC2 * sigmaValues.dsigmCa21dx * (1 + dStdDe2);
                        forceY += dSplineDCNC2 * sigmaValues.dsigmCa21dy * (1 + dStdDe2);
                        forceZ += dSplineDCNC2 * sigmaValues.dsigmCa21dz * (1 + dStdDe2);
                        forceXYZ1[0] += forceX;
                        forceXYZ1[1] += forceY;
                        forceXYZ1[2] += forceZ;
                        forceXYZ2[0] -= forceX;
                        forceXYZ2[1] -= forceY;
                        forceXYZ2[2] -= forceZ;
                        break;
                    case POLAR_HYDROPHOBIC:
                        forceX = dSplineDCNC1 * sigmaValues.dsigmCa12dx * (1 + dStdDe1);
                        forceY = dSplineDCNC1 * sigmaValues.dsigmCa12dy * (1 + dStdDe1);
                        forceZ = dSplineDCNC1 * sigmaValues.dsigmCa12dz * (1 + dStdDe1);
                        forceX += dSplineDHBC2 * sigmaValues.dsigmCa21dx * (1 + dStdDe2);
                        forceY += dSplineDHBC2 * sigmaValues.dsigmCa21dy * (1 + dStdDe2);
                        forceZ += dSplineDHBC2 * sigmaValues.dsigmCa21dz * (1 + dStdDe2);
                        forceXYZ1[0] += forceX;
                        forceXYZ1[1] += forceY;
                        forceXYZ1[2] += forceZ;
                        forceXYZ2[0] -= forceX;
                        forceXYZ2[1] -= forceY;
                        forceXYZ2[2] -= forceZ;
                        break;
                    case TWO_HYDROPHOBIC:
                        forceX  = dSplineDCNC1 * sigmaValues.dsigmCa12dx * (1 + dStdDe1);
                        forceY  = dSplineDCNC1 * sigmaValues.dsigmCa12dy * (1 + dStdDe1);
                        forceZ  = dSplineDCNC1 * sigmaValues.dsigmCa12dz * (1 + dStdDe1);
                        forceX += dSplineDCNC2 * sigmaValues.dsigmCa21dx * (1 + dStdDe2);
                        forceY += dSplineDCNC2 * sigmaValues.dsigmCa21dy * (1 + dStdDe2);
                        forceZ += dSplineDCNC2 * sigmaValues.dsigmCa21dz * (1 + dStdDe2);
                        forceXYZ1[0] += forceX;
                        forceXYZ1[1] += forceY;
                        forceXYZ1[2] += forceZ;
                        forceXYZ2[0] -= forceX;
                        forceXYZ2[1] -= forceY;
                        forceXYZ2[2] -= forceZ;
                        break;
                    case TWO_POLARS:
                        sigmaValuesBP = (SolvationDistanceAttributePolars) sigmaValues;
                        if (sigmaValuesBP.sigmHB > 0.0) { // The HB related derivatives
                            sigmaValuesBP.hydrogenBond.applyForcesToAtoms(
                                    dSplineDHBC1 * sigmaValuesBP.saltBridgeFactorA1 * (1 + dStdDe1) +
                                            dSplineDHBC2 * sigmaValuesBP.saltBridgeFactorA2  * (1 + dStdDe2));
                        }
                        break;
                }
            }
        }
    }   // Of second iteration

    public void assignForcesToAtoms() {
        int cc;
        Atom atom;
        // Finally, the appropriate forces are assigned for every atom.
        for (cc = 0; cc < molecularSystem.size(); cc++) {
            if (!molecularSystem.get(cc).atom.nowhere()) {
                atom = molecularSystem.get(cc).atom;
                if (!atom.frozen()) {
                    atom.addToFx(-forceXYZ[cc][0]); // Negating so that it is realy force (and not a mere derivative)
                    atom.addToFy(-forceXYZ[cc][1]); // Negating so that it is realy force
                    atom.addToFz(-forceXYZ[cc][2]); // Negating so that it is realy force
                }
            }
        }
    }

    public void test(TotalEnergy totalEnergy, Atom atom)  {
        double DX = 1e-7;
        double eSave,e;
        debug = false;

        System.out.println("testing Solvate");

        double cnc1 = CNC[atom.number()];
        double hbc1 = HBC[atom.number()];
        Spline2D spline1 = splines[atom.number()];
        spline1.calc(cnc1,hbc1);
        System.out.println("CNC = " + cnc1);
        System.out.println("HBC = " + hbc1);
        System.out.println("spline = " + spline1);
        double cnc23 = CNC[55];
        double hbc23 = HBC[55];
        Spline2D spline23 = splines[55];
        spline23.calc(cnc23,hbc23);
        System.out.println("CNC23 = " + cnc23);
        System.out.println("HBC23 = " + hbc23);
        System.out.println("spline23 = " + spline23);
        double s1 = spline1.s;
        double s23 = spline23.s;
        atom.addToX(DX);
            totalEnergy.evaluate();
        double cnc11 = CNC[atom.number()];
        double hbc11 = HBC[atom.number()];
        spline1 = splines[atom.number()];
        spline1.calc(cnc11,hbc11);
        System.out.println("CNC = " + cnc11);
        System.out.println("HBC = " + hbc11);
        System.out.println("spline = " + spline1);
        double cnc231 = CNC[55];
        double hbc231 = HBC[55];
        spline23 = splines[55];
        spline23.calc(cnc231,hbc231);
        System.out.println("CNC23 = " + cnc231);
        System.out.println("HBC23 = " + hbc231);
        System.out.println("spline23 = " + spline23);
        atom.addToX(-DX);
        System.out.println("Numerical splines " + (spline1.s-s1)/(cnc11-cnc1)+" "+(spline23.s-s23)/(cnc231-cnc23)+" "+(cnc231-cnc23));
        System.out.println(spline23.test(cnc23,hbc23));
        System.out.println(spline23.test(cnc231,hbc231));
        System.out.println("Neighbors:");
        for (DistanceList dl : totalEnergy.distanceMatrix().nonBondedList())
            for (Distance distance : dl)
                if ((distance.atom1() == atom) || (distance.atom2() == atom))
                    System.out.println(distance+" "+distance.atom1().type()+" "+distance.atom2().type());
        System.out.println("End of Neighbors:");


        double max = -1000;
        for (Atom atom2 : atomList) {
            double cnc = CNC[atom2.number()];
            double hbc = HBC[atom2.number()];
            Spline2D spline = splines[atom2.number()];
            spline.calc(cnc,hbc);
            if (spline.s > max) {
                max = spline.s;
            }
            String testString = spline.test(cnc,hbc);
        }

        energy = 0;
            totalEnergy.update();
        resetAuxilaryArrays();
        firstPassOverTheNonBondedList(true);
        calculateAtomEnergy(atom.core,false);
        double e1 = energy;
        energy = 0;
        atom.addToX(DX);
            totalEnergy.update();
        resetAuxilaryArrays();
        firstPassOverTheNonBondedList(true);
        double cnc2 = CNC[atom.number()];
        double hbc2 = HBC[atom.number()];
        calculateAtomEnergy(atom.core,false);
        double e2 = energy;
        double numerical = (e2-e1)/DX;
        Utils.printDebug(this, "Testing atom "+atom+" numerical = "+numerical+" cnc1 = "+cnc1+" cnc2 = "+cnc2+" hbc1 = "+hbc1+" hbc2 = "+hbc2);
        super.test(totalEnergy, atom);


    }

    private static class PerAtomWeights {
        public final double weight, sideChainCarbonWeight, polarSideChainCarbonWeight, sideChainPolarWeight, polarBackboneWeight ;
        public  PerAtomWeights(double weight, double sideChainCarbonWeight, double sideChainPolarWeight, double polarBackboneWeight, double polarSideChainCarbonWeight) {
            this.weight = weight;
            this.sideChainCarbonWeight = sideChainCarbonWeight;
            this.sideChainPolarWeight = sideChainPolarWeight;
            this.polarBackboneWeight = polarBackboneWeight;
            this.polarSideChainCarbonWeight = polarSideChainCarbonWeight;
        }
    }

}


/*
  This is what I put in the spline parameters:

  cp scModeling/newForClust/SolvateCarbonSplinesCorrection1.txt ../meshi/parameters/meshiPotential/SolvateCarbonSideChainSplines.dat

  cp scModeling/newForClust/SolvatePolarBackboneSplinesOldFashion_noCorr.txt ../meshi/parameters/meshiPotential/SolvatePolarBackboneSplines.dat

  cp scModeling/newForClust/SolvatePolarBackboneSplines_NOTBAYES_NoHB.txt ../meshi/parameters/meshiPotential/SolvatePolarBackboneSplines_alt2.dat

  cp scModeling/newForClust/SolvatePolarSideChainSplinesOldFashion_noCorr.txt ../meshi/parameters/meshiPotential/SolvatePolarSideChainSplines.dat

  cp scModeling/newForClust/SolvatePolarSideChainSplines_NOTBAYES_NoHB.txt ../meshi/parameters/meshiPotential/SolvatePolarSideChainSplines_alt2.dat

*/
