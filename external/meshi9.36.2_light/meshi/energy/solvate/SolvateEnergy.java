/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.energy.CooperativeEnergyTerm;
import meshi.energy.EnergyInfoElement;
import meshi.energy.TotalEnergy;
import meshi.energy.solvate.hydrogenBonds.AbstractHydrogenBondList;
import meshi.energy.solvate.hydrogenBonds.HydrogenBondDahiyatList;
import meshi.geometry.Distance;
import meshi.geometry.DistanceList;
import meshi.geometry.DistanceLists;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;
import meshi.util.Stat;
import meshi.util.Utils;
import meshi.util.mathTools.Spline1D;
import meshi.util.mathTools.Spline2D;

/**
 * The implementation of the cooperative solvation term for proteins as described in Kalisman & Keasar (2008).
 * Since the cooperative solvation is described in the above paper, we bring here only the implementaion
 * details. Especially regarding the calculation of the derivatives, which was too lengthy for the paper.
 * <p/>
 * The class allows to give a different weight to the solvation energies of certain atom types in the final
 * summation described in Eq. 5. These atom types are side-chain polars, side-chain carbons, and backbone
 * polars. The class also include a regular hydrogen bond term. The functional form is therefore:
 * <p/>
 * Esolv = weightSCPolarSolvate*Eside_chain_polars +
 * weightSCCarbonSolvate*Eside_chain_carbons +
 * weightBBPolarSolvate*Ebackbone_polars +
 * weightHB*Ehb
 * <p/>
 * Where Ehb is the negative of the HBC summation over all atoms. The weights are defined in the "Creator" class,
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
 * with respect to the distance atoms are alsoclaculated.  Since some of this values will be needed also in step 3,
 * we save them in an instance of "SolvateDistanceAttribute" that is attached as an "Attribute" to the Distance instance.
 * 2) A pass isOn the atom list. Once we have the CNC and HBC of every atom in the protein, we can proceed to calculate the solvation energy
 * associated with every atom. In this implementation we combined the EI(CNC,HBC) evaluation (Eq. 3) of every atom and
 * the -log(spline(EI)) evaluation (Eq. 4) into a single step by using a 2D spline, i.e. spline2D(CNC,HBC). The 2D spline
 * is, of course, atom type specific. The derivatives of each atom solvate energy value with respect to the HBC and CNC
 * are also calculated.
 * 3) A second pass over the non-bonded list. Equiped with the The derivatives of the atom energies (with respect
 * to the CNC's and HBC's) from step 2, we can now calculate the energy derivative with respect to the atomic
 * coordinates. In this step we simply make sure that every term that arises from the derivative chain rule is accounted for.
 */
public final class SolvateEnergy extends CooperativeEnergyTerm {

    // Relative strength of SALT BRIDGES compared with regular HYDROGEN BONDS for desolvation purposes, i.e. in
    // regard to the effect isOn observed CNC medians. Following Table 1 in the paper.
    public static final double SALT_BRIDGE_STRENGTH_ASP_OD = 0.45;
    public static final double SALT_BRIDGE_STRENGTH_GLU_OE = 0.45;
    public static final double SALT_BRIDGE_STRENGTH_LYS_NZ = 1.65;
    public static final double SALT_BRIDGE_STRENGTH_ARG_NH = 1.5;
    public static final double SALT_BRIDGE_STRENGTH_TRO = 0.45;
    public static final double SALT_BRIDGE_STRENGTH_TRN = 1.8;
    // A flag for distances that do not take part in this term
    public static final SolvateDistanceAttribute NOT_PATRICIPATING_DISTANCE = new SolvateDistanceAttribute();
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
    private double[] HBC;  // Eq. 2
    public double[] dSplineDCNC;
    public double[] dSplineDHBC;
    public double[] forceX;
    public double[] forceY;
    public double[] forceZ;

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
        CARBON, BACKBONE_POLAR, SIDECHAIN_POLAR, HYDROGEN
    }

    private SuperType[] superType;
    /**
     * The 2D spline array (i.e. spline(CNC,HBC)) for the polar side-chain atoms. The indexing in the array is the
     * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) *
     */
    private Spline2D[] splinesSCPolar;
    /**
     * The 2D spline array (i.e. spline(CNC,HBC)) for the polar backbone atoms. The indexing in the array is the
     * index of the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded) *
     */
    private Spline2D[] splinesBB;
    /**
     * The 1D spline array (i.e. spline(CNC,HBC)) for the polar carbon atoms. The indexing in the array is the index of
     * the atom type (i.e. 0-189, because we have 190 atom types in meshi currently. The number 190 is not hard coded). The splines
     * are 1D because HBC is 0, for carbons. *
     */
    private Spline1D[] splinesSCCarbon;
    /**
     * The only hydrogen bond list class in the term. *
     */
    public AbstractHydrogenBondList solvateHB;


    // Weights
    private double weightSCPolarSolvate;
    private double weightBBPolarSolvate;
    private double weightSCCarbonSolvate;

    public final MolecularSystem molecularSystem;
    public final int molecularSystemSize;
    private Atom atom;
    private int atomTypeNumber;
    private AtomType atomType;
    private Stat stat;

    public SolvateEnergy() {
        molecularSystem = null;
        molecularSystemSize = -1;
    }

    public double[] atomEnergies;

    /**
     * See the comment at the top of the class for descriptions isOn the weights.
     */

    public SolvateEnergy(AtomList atomList,
                         DistanceMatrix dm,
                         SolvateParametersList parameters,
                         EnergyInfoElement info) {
        super(toArray(), atomList, dm, parameters, info);
        SolvateInfoElement solvateInfo = (SolvateInfoElement) info;

        molecularSystem = atomList.molecularSystem();

        weightSCPolarSolvate  = ((EnergyInfoElement) solvateInfo.scPolarInfo).weight();
        weightBBPolarSolvate  = ((EnergyInfoElement) solvateInfo.bbPolarInfo).weight();
        weightSCCarbonSolvate = ((EnergyInfoElement) solvateInfo.scCarbonInfo).weight();

        this.comment = "Solvation";
        int c;

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
        dSplineDCNC = new double[molecularSystemSize];
        dSplineDHBC = new double[molecularSystemSize];
        forceX = new double[molecularSystemSize];
        forceY = new double[molecularSystemSize];
        forceZ = new double[molecularSystemSize];
        superType = new SuperType[molecularSystemSize];

        // Setting up the splines array, that can associate each residue in the protein with the spline that
        // suits its type.
        splinesSCPolar = new Spline2D[molecularSystemSize];
        splinesSCCarbon = new Spline1D[molecularSystemSize];
        splinesBB = new Spline2D[molecularSystemSize];
        for (c = 0; c < molecularSystemSize; c++) {
            atom = molecularSystem.get(c).atom;
            if (!atom.nowhere()) {
                atomTypeNumber = atom.type().ordinal();
                splinesSCPolar[c] = parameters.scPolarSplines[atomTypeNumber];
                splinesSCCarbon[c] = parameters.scCarbonSplines[atomTypeNumber];
                splinesBB[c] = parameters.bbSplines[atomTypeNumber];
            }
        }

        // Setting up the HB list
        solvateHB = new HydrogenBondDahiyatList(dm, atomList, parameters);
        updateableResources.add(solvateHB);

        // Determining the superType for each atom: (0) Hydrophobic, (1) Backbone polar, (2) Sidechain polar (3) Hydrogens
        for (c = 0; c < molecularSystemSize; c++) {
            atom = molecularSystem.get(c).atom;
            atomType = atom.type();
            if (atomType.isCarbon() || (atomType == AtomType.MSD))
                superType[c] = SuperType.CARBON;
            else if ((atomType.isOxygen() || atomType.isNitrogen()) && (atom.name().length() == 1))
                superType[c] = SuperType.BACKBONE_POLAR;
            else if ((atomType.isOxygen() || atomType.isNitrogen() || (atomType == AtomType.CSG))
                    && (atom.name().length() > 1))
                superType[c] = SuperType.SIDECHAIN_POLAR;
            else if (atomType.isHydrogen())
                superType[c] = SuperType.HYDROGEN;
            else
                throw new RuntimeException("The following if's should have covered all cases. But due to my carelessness this atom was left out:\n" +
                        atom + "\n");
            stat = new Stat();
        }

    }  // Of constructor

    public void evaluateAtoms() {
        evaluate(true, weightSCPolarSolvate, weightSCCarbonSolvate, weightBBPolarSolvate);
    }

    /**
     * Calculates Esolv with the weights given in the constructor.
     */
    public EnergyInfoElement evaluate() {
        return evaluate(false, weightSCPolarSolvate, weightSCCarbonSolvate, weightBBPolarSolvate);
    }

    /**
     * Calculates Esolv with the weights you give as parameters!
     */
    public final EnergyInfoElement evaluate(boolean updateAtoms, double W_SCPolarSolvate, double W_SCCarbonSolvate, double W_BBPolarSolvate) {

        double energy;
        stat.reset();
        resetAuxilaryArrays();
        firstPassOverTheNonBondedList();
        energy = calculatingTheSolvationEnergiesOfEachAtom(updateAtoms, W_SCPolarSolvate, W_SCCarbonSolvate, W_BBPolarSolvate);
        secondPassOverTheNonBondedList();
        assignForcesToAtoms();
        ((SolvateInfoElement) info()).stdElememt().setValue(stat.getStd());
        info.setValue(energy);
        return info;
    }

    private void resetAuxilaryArrays() {
        int cc;
        //Reseting the auxilary arrays
        for (cc = 0; cc < molecularSystem.size(); cc++) {
            CNC[cc] = 0;
            HBC[cc] = 0;
            dSplineDCNC[cc] = 0;
            dSplineDHBC[cc] = 0;
            forceX[cc] = forceY[cc] = forceZ[cc] = 0.0;
        }
    }

    public void firstPassOverTheNonBondedList() {
        int atomNumber1, atomNumber2;
        SolvateDistanceType type;

        DistanceLists dislist = dm.nonBondedList();
        SolvateDistanceAttribute sigmaValues = null;
        SolvateDistanceAttributeBetweenPolars sigmaValuesBP = null;
        SolvateDistanceAttributeWithNonPolar sigmaValuesWNP = null;
        SolvateDistanceAttributeNonPolarPolar sigmaValuesNPP = null;
        SolvateDistanceAttributePolarNonPolar sigmaValuesPNP = null;
        // ***********************************
        // First pass over the non-bonded list
        // ***********************************
        for (DistanceList distanceList : dislist) {
            for (Distance dis : distanceList) {
                sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
                // In case there is not a "SolvateDistanceAttribute" in the distance, we create one.
                if (sigmaValues == null) {
                    sigmaValues = SolvateDistanceAttributeCreator.create(dis, parameters, solvateHB);
                    dis.addAttribute(sigmaValues);
                }
                if (sigmaValues == NOT_PATRICIPATING_DISTANCE) continue;
                type = sigmaValues.type();
                //Distances with hydrogens are ignored.
                atomNumber1 = dis.atom1.number;
                atomNumber2 = dis.atom2.number;
                switch (type) {
                    case HYDROPHOBIC_POLAR:
                        sigmaValuesNPP = (SolvateDistanceAttributeNonPolarPolar) sigmaValues;
                        sigmaValuesNPP.update();
                        CNC[atomNumber2] += sigmaValuesNPP.sigmCa2;
                        break;
                    case POLAR_HYOPHOBIC:
                        sigmaValuesPNP = (SolvateDistanceAttributePolarNonPolar) sigmaValues;
                        sigmaValuesPNP.update();
                        CNC[atomNumber1] += sigmaValuesPNP.sigmCa1;
                        break;
                    case TWO_HYDROPHOBIC:
                        sigmaValuesWNP = (SolvateDistanceAttributeWithNonPolar) sigmaValues;
                        sigmaValuesWNP.update();
                        CNC[atomNumber1] += sigmaValuesWNP.sigmCa1;
                        CNC[atomNumber2] += sigmaValuesWNP.sigmCa2;
                        break;
                    case TWO_POLARS:
                        sigmaValuesBP = (SolvateDistanceAttributeBetweenPolars) sigmaValues;
                        sigmaValuesBP.update();
                        HBC[atomNumber1] += sigmaValuesBP.saltBridgeFactorA1 * sigmaValuesBP.sigmHB;
                        HBC[atomNumber2] += sigmaValuesBP.saltBridgeFactorA2 * sigmaValuesBP.sigmHB;
                        break;
                }
            }
        }
    }


    public double calculatingTheSolvationEnergiesOfEachAtom(boolean updateAtoms, double W_SCPolarSolvate, double W_SCCarbonSolvate, double W_BBPolarSolvate) {
        int atomNumber;
        double energy = 0;
        double energySCPolar;
        double energySCCarbon;
        double energyBBON;
        double energySCPolarTotal = 0;
        double energySCCarbonTotal = 0;
        double energyBBONtotal = 0;
        // ***************************************************************************
        // A pass over the atom list. Calculating the solvation energies of each atom.
        // ***************************************************************************
        if (W_SCCarbonSolvate > 300)
            throw new RuntimeException(""+W_SCCarbonSolvate);
        for (int i = 0; i < atomEnergies.length; i++) atomEnergies[i] = 0;

        for (AtomCore atom : molecularSystem) {
            if (atom.atom.nowhere()) continue;
            energySCPolar = energySCCarbon = energyBBON = 0.0;
            atomNumber = atom.number;
            switch (superType[atomNumber]) {
                case HYDROGEN:
                    break;
                case CARBON:
                    try {
                        splinesSCCarbon[atomNumber].calc(CNC[atomNumber]);
                    } catch (Exception e) {
                        throw new RuntimeException("SCC: " + CNC[atomNumber] + " " + HBC[atomNumber] + "\n" +
                                molecularSystem.get(atomNumber).atom + "\n" + e);
                    }
                    energySCCarbon = W_SCCarbonSolvate * splinesSCCarbon[atomNumber].s;
                    dSplineDCNC[atomNumber] = W_SCCarbonSolvate * splinesSCCarbon[atomNumber].s_tag;
                    dSplineDHBC[atomNumber] = 0.0;
                    break;
                case BACKBONE_POLAR:
                    try {
                        splinesBB[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]);
                    } catch (Exception e) {
                        throw new RuntimeException("BBON: " + CNC[atomNumber] + " " + HBC[atomNumber] + "\n" +
                                molecularSystem.get(atomNumber).atom + "\n");
                    }
                    energyBBON = W_BBPolarSolvate * (splinesBB[atomNumber].s);
                    dSplineDCNC[atomNumber] = W_BBPolarSolvate * splinesBB[atomNumber].s_tag_x;
                    dSplineDHBC[atomNumber] = W_BBPolarSolvate * splinesBB[atomNumber].s_tag_y;
                    break;
                case SIDECHAIN_POLAR:
                    try {
                        splinesSCPolar[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]);
                    } catch (Exception e) {
                        throw new RuntimeException("SCP: " + CNC[atomNumber] + " " + HBC[atomNumber] + "\n" +
                                molecularSystem.get(atomNumber).atom + "\n");
                    }
                    energySCPolar = W_SCPolarSolvate * splinesSCPolar[atomNumber].s;
                    dSplineDCNC[atomNumber] = W_SCPolarSolvate * splinesSCPolar[atomNumber].s_tag_x;
                    dSplineDHBC[atomNumber] = W_SCPolarSolvate * splinesSCPolar[atomNumber].s_tag_y;
                    break;
            }

            energySCPolarTotal += energySCPolar;
            energySCCarbonTotal += energySCCarbon;
            energyBBONtotal += energyBBON;
            stat.add(energySCPolar + energySCCarbon + energyBBON);
            energy += (energySCPolar + energySCCarbon + energyBBON); // We want to count each HB once.
            if (updateAtoms)
                molecularSystem.get(atomNumber).atom.addEnergy(energySCPolar + energySCCarbon + energyBBON);
            if (molecularSystem.get(atomNumber).atom.number() != atomNumber)
                throw new RuntimeException("Check the size of atomEnergies");
            atomEnergies[molecularSystem.get(atomNumber).atom.number()] += energySCPolar + energySCCarbon + energyBBON;
        }
        ( info.getChildren().get(0)).setValue(energySCPolarTotal);
        ( info.getChildren().get(1)).setValue(energySCCarbonTotal);
        ( info.getChildren().get(2)).setValue(energyBBONtotal);
         return energy;
    }

    public void secondPassOverTheNonBondedList() {
        int atomNumber1, atomNumber2;
        SolvateDistanceType type;
        DistanceLists dislist = dm.nonBondedList();
        SolvateDistanceAttribute sigmaValues = null;
        SolvateDistanceAttributeBetweenPolars sigmaValuesBP = null;
        SolvateDistanceAttributeWithNonPolar sigmaValuesWNP = null;
        SolvateDistanceAttributeNonPolarPolar sigmaValuesNPP = null;
        SolvateDistanceAttributePolarNonPolar sigmaValuesPNP = null;
        // ************************************
        // Second pass over the non-bonded list
        // ************************************
        double sumForce;
        double force;
        for (DistanceList distanceList : dislist) {
            for (Distance dis : distanceList) {
                sigmaValues = (SolvateDistanceAttribute) dis.getAttribute(SolvateDistanceAttribute.SOLVATE_ALL_ATOM_ATTRIBUTE);
                if (sigmaValues == NOT_PATRICIPATING_DISTANCE) continue;
                type = sigmaValues.type();
                atomNumber1 = dis.atom1.number;
                atomNumber2 = dis.atom2.number;
                switch (type) {
                    case HYDROPHOBIC_POLAR:
                        sumForce = 0;
                        sigmaValuesNPP = (SolvateDistanceAttributeNonPolarPolar) sigmaValues;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dx;
                        forceX[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dy;
                        forceY[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dz;
                        forceZ[atomNumber1] += force;
                        sumForce += force;
                        forceX[atomNumber2] -= dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dx;
                        forceY[atomNumber2] -= dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dy;
                        forceZ[atomNumber2] -= dSplineDCNC[atomNumber2] * sigmaValuesNPP.dsigmCa2dz;
                        //if (Math.abs(sumForce) > 50) System.out.println("\nHot solvate contact HYDROPHOBIC_POLAR\n"+dis+"\n"+"Force = "+sumForce);
                        break;
                    case POLAR_HYOPHOBIC:
                        sumForce = 0;
                        sigmaValuesPNP = (SolvateDistanceAttributePolarNonPolar) sigmaValues;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dx;
                        forceX[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dy;
                        forceY[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dz;
                        forceZ[atomNumber1] += force;
                        sumForce += force;
                        forceX[atomNumber2] -= dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dx;
                        forceY[atomNumber2] -= dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dy;
                        forceZ[atomNumber2] -= dSplineDCNC[atomNumber1] * sigmaValuesPNP.dsigmCa1dz;
                        // if (Math.abs(sumForce) > 50) System.out.println("\nHot solvate contact POLAR_HYOPHOBIC\n"+dis+"\n"+"Force = "+sumForce);
                        break;
                    case TWO_HYDROPHOBIC:
                        sumForce = 0;
                        sigmaValuesWNP = (SolvateDistanceAttributeWithNonPolar) sigmaValues;
                        // Doing the self derivatives
                        force = dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dx;
                        forceX[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dy;
                        forceY[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dz;
                        forceZ[atomNumber1] += force;
                        sumForce += force;
                        // if (Math.abs(sumForce) > 50) System.out.println("\nHot solvate contact TWO_HYDROPHOBIC self\n"+dis+"\n"+"Force = "+sumForce);

                        sumForce = 0;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dx;
                        forceX[atomNumber2] -= force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dy;
                        forceY[atomNumber2] -= force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dz;
                        forceZ[atomNumber2] -= force;
                        sumForce += force;
                        // if (Math.abs(sumForce) > 50) System.out.println("\nHot solvate contact TWO_HYDROPHOBIC self 2\n"+dis+"\n"+"Force = "+sumForce);

                        // Doing the cross derivatives
                        sumForce = 0;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dx;
                        forceX[atomNumber2] -= force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dy;
                        forceY[atomNumber2] -= force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber1] * sigmaValuesWNP.dsigmCa1dz;
                        forceZ[atomNumber2] -= force;
                        sumForce += force;
                        // if (Math.abs(sumForce) > 50) System.out.println("\nHot solvate contact TWO_HYDROPHOBIC cross 1\n"+dis+"\n"+"Force = "+sumForce);

                        sumForce = 0;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dx;
                        forceX[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dy;
                        forceY[atomNumber1] += force;
                        sumForce += force;
                        force = dSplineDCNC[atomNumber2] * sigmaValuesWNP.dsigmCa2dz;
                        forceZ[atomNumber1] += force;
                        sumForce += force;
                        //     if (Math.abs(sumForce) > 50) System.out.println("\nHot solvate contact TWO_HYDROPHOBIC cross 2\n"+dis+"\n"+"Force = "+sumForce);

                        break;
                    case TWO_POLARS:
                        sigmaValuesBP = (SolvateDistanceAttributeBetweenPolars) sigmaValues;
                        if (sigmaValuesBP.sigmHB > 0.0) { // The HB related derivatives
                            solvateHB.findBondByPolars(dis.atom1(), dis.atom2()).applyForcesToAtoms(
                                    dSplineDHBC[atomNumber1] * sigmaValuesBP.saltBridgeFactorA1 +
                                            dSplineDHBC[atomNumber2] * sigmaValuesBP.saltBridgeFactorA2);
                            //                              - 0.5 * W_HB * (sigmaValuesBP.saltBridgeFactorForHBenergyA1 + sigmaValuesBP.saltBridgeFactorForHBenergyA2));
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
                    atom.addToFx(-forceX[cc]); // Negating so that it is realy force (and not a mere derivative)
                    atom.addToFy(-forceY[cc]); // Negating so that it is realy force
                    atom.addToFz(-forceZ[cc]); // Negating so that it is realy force
                }
            }
        }
    }

    public void test(TotalEnergy totalEnergy, Atom atom) {
        double DX = 1e-7;
        double eSave,e;

        System.out.println("Starting test of "+this);

        double saveWeightSCPolarSolvate = weightSCPolarSolvate;
        double saveWeightSCCarbonSolvate = weightSCCarbonSolvate;
        double saveWeightBBPolarSolvate = weightBBPolarSolvate;
        System.out.println("testing SCPolarSolvate");

        weightSCCarbonSolvate = 0;
        weightBBPolarSolvate = 0;
        super.test(totalEnergy,atom);

        System.out.println("testing SCCarbonSolvate");
        weightSCPolarSolvate = 0;
        weightSCCarbonSolvate = saveWeightSCCarbonSolvate;
        super.test(totalEnergy,atom);
        double[] cncSave = new double[CNC.length];
        double[] splineSave = new double[splinesSCCarbon.length];
        double sum,sumSave;
        try {
            totalEnergy.update();
        }
        catch (Exception ex) {throw new RuntimeException(ex);}
        sumSave = 0;
        resetAuxilaryArrays();
        firstPassOverTheNonBondedList();
        eSave = calculatingTheSolvationEnergiesOfEachAtom(false, 0, 1, 0);
        for (int i = 0; i < CNC.length; i++) {
            cncSave[i] = CNC[i];
            splinesSCCarbon[i].calc(CNC[i]);
            splineSave[i] = splinesSCCarbon[i].s;
            sumSave += splineSave[i];
        }

        double[][] coordinates = new double[3][];
        coordinates[0] = atom.X();
        coordinates[1] = atom.Y();
        coordinates[2] = atom.Z();
        for (int iCoor = 0; iCoor < 3; iCoor++) {
            for (int jDirection = -1; jDirection < 3; jDirection += 2) {
                coordinates[iCoor][0] += jDirection*DX;
                sum = 0;
                try {totalEnergy.update();}
                catch (Exception ex) {throw new RuntimeException(ex);}
                resetAuxilaryArrays();
                firstPassOverTheNonBondedList();
                e = calculatingTheSolvationEnergiesOfEachAtom(false, 0, 1, 0);
                for (int i = 0; i < CNC.length; i++) {
                    if (Math.abs(cncSave[i]-CNC[i])> 0.01)  {
                        System.out.println("problematic atom "+molecularSystem.get(i));
                        System.out.println("Problem with CCarbonSolvate cnc "+iCoor+" "+jDirection+" CNC["+i+"] = "+
                                                                      CNC[i]+" cncSave = "+cncSave[i]);
                    }
                    splinesSCCarbon[i].calc(CNC[i]);
                    if (Math.abs(splineSave[i]-splinesSCCarbon[i].s)> 0.01) {
                        System.out.println("problematic atom "+molecularSystem.get(i));
                        System.out.println("Problem with CCarbonSolvate splines "+iCoor+" "+jDirection+" CNC["+i+"] = "+
                                                                      CNC[i]+" cncSave = "+cncSave[i]+"\n"+
                                                                      "splineSave = "+splineSave[i]+"\n"+splinesSCCarbon[i]);
                    }
                    sum += splinesSCCarbon[i].s;
                }
                System.out.println(iCoor+" "+jDirection+" sum1 = "+sum+"; sumSave = "+sumSave+"; e = "+e+
                                                        "; eSave "+eSave+"; (e-eSave)/DX = "+(e-eSave)/(jDirection*DX));
                coordinates[iCoor][0] -= jDirection*DX;
            }
        }


        System.out.println("testing BBPolarSolvate");
        weightSCCarbonSolvate = 0;
        weightBBPolarSolvate = saveWeightBBPolarSolvate;
        super.test(totalEnergy,atom);

        weightSCCarbonSolvate = saveWeightSCCarbonSolvate;
        weightSCPolarSolvate  = saveWeightSCPolarSolvate;

        try {totalEnergy.update();}
        catch (Exception ex) {throw new RuntimeException(ex);}
        resetAuxilaryArrays();
        firstPassOverTheNonBondedList();
        for (AtomCore atom1 : molecularSystem) {
            if (atom1.atom.nowhere()) continue;
            int atomNumber = atom1.number;
            splinesSCCarbon[atomNumber].calc(CNC[atomNumber]);
            double s1 = splinesSCCarbon[atomNumber].s;
            splinesSCCarbon[atomNumber].calc(CNC[atomNumber]+DX);
            double s2 = splinesSCCarbon[atomNumber].s;
            double s3;
            if(CNC[atomNumber] > DX){
                splinesSCCarbon[atomNumber].calc(CNC[atomNumber]-DX);
                s3 = splinesSCCarbon[atomNumber].s;
            } else s3 = s1;
            if ((Math.abs(splinesSCCarbon[atomNumber].s_tag-(s2-s1)/DX)>0.1) ||
                 (Math.abs(splinesSCCarbon[atomNumber].s_tag-(s1-s3)/DX)>0.1)) {
                System.out.println("Problem with splinesSCCarbon "+atom1);
                System.out.println("CNC = "+ CNC[atomNumber]);
                System.out.println("s1 = "+s1+"; s2 = "+s2+"; s3 = "+s3+"; s_tag = "+
                                   splinesSCCarbon[atomNumber].s_tag+"; (s2-s1)/DX = "+(s2-s1)/DX);
                System.out.println(splinesSCCarbon[atomNumber]);
            }

            splinesBB[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]);
            s1 = splinesBB[atomNumber].s;
            splinesBB[atomNumber].calc(CNC[atomNumber]+DX, HBC[atomNumber]);
             s2 = splinesBB[atomNumber].s;
            if (Math.abs(splinesBB[atomNumber].s_tag_x-(s2-s1)/DX)>0.1) {
                System.out.println("Problem with splinesBB x "+atom1);
                System.out.println("CNC = "+ CNC[atomNumber]);
                System.out.println("s1 = "+s1+"; s2 = "+s2+"; s_tag_x = "+
                                   splinesBB[atomNumber].s_tag_x+"; (s2-s1)/DX = "+(s2-s1)/DX);
            }

            splinesBB[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]);
                    s1 = splinesBB[atomNumber].s;
                    splinesBB[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]+DX);
                     s2 = splinesBB[atomNumber].s;
                    if (Math.abs(splinesBB[atomNumber].s_tag_y-(s2-s1)/DX)>0.1) {
                        System.out.println("Problem with splinesBB y "+atom1);
                        System.out.println("CNC = "+ CNC[atomNumber]);
                        System.out.println("s1 = "+s1+"; s2 = "+s2+"; s_tag_y = "+
                                           splinesBB[atomNumber].s_tag_y+"; (s2-s1)/DX = "+(s2-s1)/DX);
                    }

            splinesSCPolar[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]);
            s1 = splinesSCPolar[atomNumber].s;
            splinesSCPolar[atomNumber].calc(CNC[atomNumber]+DX, HBC[atomNumber]);
             s2 = splinesSCPolar[atomNumber].s;
            if (Math.abs(splinesSCPolar[atomNumber].s_tag_x-(s2-s1)/DX)>0.1) {
                System.out.println("Problem with splinesSCPolar x "+atom1);
                System.out.println("CNC = "+ CNC[atomNumber]);
                System.out.println("s1 = "+s1+"; s2 = "+s2+"; s_tag_x = "+
                                   splinesSCPolar[atomNumber].s_tag_x+"; (s2-s1)/DX = "+(s2-s1)/DX);
            }

            splinesSCPolar[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]);
            s1 = splinesSCPolar[atomNumber].s;
            splinesSCPolar[atomNumber].calc(CNC[atomNumber], HBC[atomNumber]+DX);
             s2 = splinesSCPolar[atomNumber].s;
            if (Math.abs(splinesSCPolar[atomNumber].s_tag_y-(s2-s1)/DX)>0.1) {
                System.out.println("Problem with splinesSCPolar y "+atom1);
                System.out.println("CNC = "+ CNC[atomNumber]);
                System.out.println("s1 = "+s1+"; s2 = "+s2+"; s_tag_y = "+
                                   splinesSCPolar[atomNumber].s_tag_y+"; (s2-s1)/DX = "+(s2-s1)/DX);
            }
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
