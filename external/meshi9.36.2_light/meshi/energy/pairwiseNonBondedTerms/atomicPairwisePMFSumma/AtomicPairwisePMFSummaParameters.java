/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma;

import meshi.parameters.AtomType;
import meshi.parameters.MeshiPotential;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.file.MeshiLineReader;

import java.util.Arrays;
import java.util.StringTokenizer;

/**
 * Parameters for the atomic pairwise PMF. Class includes polynomial splines
 * for all atom-atom PMFs.
 *
 * @author El-ad David Amir
 */
public class AtomicPairwisePMFSummaParameters implements KeyWords, MeshiPotential {

    /* labels used as part of the parameters file format */
    private static final String LABEL_BINS_BEGIN = "BINS_BEGIN";
    private static final String LABEL_BINS_END = "BINS_END";
    private static final String LABEL_ATOM_NAMES_BEGIN = "ATOM_NAMES_BEGIN";
    private static final String LABEL_ATOM_NAMES_END = "ATOM_NAMES_END";
    private static final String LABEL_PMF_BEGIN = "PMF_BEGIN";
    private static final String LABEL_PMF_END = "PMF_END";
    private static final int[] MESHI2SUMMA = new int[AtomType.values().length];
    /* internal variable guaranteeing parameters are only loaded once */
    private boolean parametersLoaded = false;
    /* list of bins of spline */
    private double[] bins;
    /* list of Summa atom names; names order correspond to PMFs order */
    private String[] atomNames;
    /* coefficients for all splines */
    private CoefficientsMatrixForAtomPairSpline[][] coefs;

    /**
     * Load parameters list from file.
     */
    public AtomicPairwisePMFSummaParameters(CommandList commands) {
        /* verify parameters are only loaded once */
        if (parametersLoaded)
            return;
        parametersLoaded = true;

        /* open parameters file for input */
        String parametersFileName = commands.firstWord(PARAMETERS_DIRECTORY).secondWord() + "/" + ATOMIC_PAIRWISE_PMF_SUMMA_PARAMETERS;
        MeshiLineReader mlrParameters = new MeshiLineReader(parametersFileName);

        String line;

        /* read bins */
        /* read LABEL_BINS_BEGIN */
        line = mlrParameters.readLine("#");
        if (!line.equals(LABEL_BINS_BEGIN))
            loadingParametersFailed(line);
        /* read nu. of bins */
        line = mlrParameters.readLine("#");
        bins = new double[Integer.valueOf(line).intValue()];
        /* read bins */
        for (int i = 0; i < bins.length; i++) {
            line = mlrParameters.readLine("#");
            bins[i] = Double.valueOf(line).doubleValue();
        }
        /* read LABEL_BINS_END */
        line = mlrParameters.readLine("#");
        if (!line.equals(LABEL_BINS_END))
            loadingParametersFailed(line);

        /* read atom names */
        /* read LABEL_ATOM_NAMES_BEGIN */
        line = mlrParameters.readLine("#");
        if (!line.equals(LABEL_ATOM_NAMES_BEGIN))
            loadingParametersFailed(line);
        /* read nu. of atom names */
        line = mlrParameters.readLine("#");
        atomNames = new String[Integer.valueOf(line).intValue()];
        /* read atom names */
        for (int i = 0; i < atomNames.length; i++) {
            line = mlrParameters.readLine("#");
            atomNames[i] = line;
        }
        /* read LABEL_ATOM_NAMES_END */
        line = mlrParameters.readLine("#");
        if (!line.equals(LABEL_ATOM_NAMES_END))
            loadingParametersFailed(line);

        /* read coefficients; a coefs matrix for each atom pair */
        coefs = new CoefficientsMatrixForAtomPairSpline[atomNames.length][atomNames.length];
        int numCoefsRead = 0;        /* number of coefficients matrices read */
        /* continue reading coefs till end of file is reached */
        while ((line = mlrParameters.readLine("#")) != null) {
            /* read LABEL_PMF_BEGIN */
            if (!line.equals(LABEL_PMF_BEGIN))
                loadingParametersFailed(line);
            /* read atom numbers; conversion from matlab format to java format,
            * we substract one from each atom number */
            int atom1Number, atom2Number;
            line = mlrParameters.readLine("#");
            atom1Number = Integer.valueOf(line).intValue() - 1;
            line = mlrParameters.readLine("#");
            atom2Number = Integer.valueOf(line).intValue() - 1;
            /* read coefs using string tokenizer; each line has four coefs, for one bin*/
            coefs[atom1Number][atom2Number] = new CoefficientsMatrixForAtomPairSpline(bins, bins.length - 1,
                    atomNames[atom1Number] + "," + atomNames[atom2Number]);
            for (int i = 0; i < bins.length - 1; i++) {
                line = mlrParameters.readLine("#");
                StringTokenizer st = new StringTokenizer(line);
                coefs[atom1Number][atom2Number].coefs[i][0] = Double.valueOf(st.nextToken()).doubleValue();
                coefs[atom1Number][atom2Number].coefs[i][1] = Double.valueOf(st.nextToken()).doubleValue();
                coefs[atom1Number][atom2Number].coefs[i][2] = Double.valueOf(st.nextToken()).doubleValue();
                coefs[atom1Number][atom2Number].coefs[i][3] = Double.valueOf(st.nextToken()).doubleValue();
            }
            /* read LABEL_PMF_END */
            line = mlrParameters.readLine("#");
            if (!line.equals(LABEL_PMF_END))
                loadingParametersFailed(line);

            numCoefsRead++;
        }
        for (AtomType type : AtomType.values())
            MESHI2SUMMA[type.ordinal()] = summaAtomNameToParametersAtomIndex(meshiAtomTypeToSummaAtomName(type));
    }

    /**
     * Throws exception. Called when loading the file failed.
     */
    private void loadingParametersFailed(String line) {
        throw new RuntimeException("unable to load atomic pairwise PMF parameters; read the line:\n" + line);
    }

    /**
     * Retrieve spline for atom pair. Return null if not found.
     */
    public CoefficientsMatrixForAtomPairSpline splineForAtomPair(AtomType atom1Type,
                                                                 AtomType atom2Type) {
        return splineForAtomPair(atom1Type.ordinal(), atom2Type.ordinal());
    }

    public CoefficientsMatrixForAtomPairSpline splineForAtomPair(int atom1TypeOrdinal, int atom2TypeOrdinal) {
        int atom1Index = MESHI2SUMMA[atom1TypeOrdinal];
        int atom2Index = MESHI2SUMMA[atom2TypeOrdinal];

        if (atom1Index == -1 || atom2Index == -1)
            return null;

        /* only half of the coefficients are stored (since the distance between a to b is
       * identical to the distance between b to a). make sure we retrieve coefficients
       * that were actually calculated.
       */
// 	if( atom1Index <= atom2Index ) {
// 	    if  (coefs[atom1Index][atom2Index] == null) throw new RuntimeException("Failed to find index1\n"+atom1Type+" "+atom1Index+"\n"+atom2Type+" "+atom2Index);
// 	}
// 	else if  (coefs[atom2Index][atom1Index] == null) throw new RuntimeException("Failed to find index2\n"+atom1Type+" "+atom1Index+"\n"+atom2Type+" "+atom2Index);

        if (atom1Index <= atom2Index) {
            return coefs[atom1Index][atom2Index];
        } else {
            return coefs[atom2Index][atom1Index];
        }
    }

    /**
     * Translates MESHI atom type to Summa atom name. Return null if not found.
     */
    private static String meshiAtomTypeToSummaAtomName(AtomType atomType) {
        int index = Arrays.binarySearch(atomNamesTranslationsMESHI, atomType.toString());

        if (index < 0)
            return null;

        return atomNamesTranslationsSumma[index];
    }

    /**
     * Translates Summa atom name to parameters file atom index. Return -1 if not found.
     */
    private int summaAtomNameToParametersAtomIndex(String name) {
        for (int index = 0; index < atomNames.length; index++)
            if (atomNames[index].equals(name))
                return index;

        return -1;
    }

    /* the Summa energy term uses different atom names than MESHI. these two
     * arrays are used in the translation from MESHI names to Summa names, and
     * vice versa. index i in both arrays corresponds to the same atom type.
     */
    private static final String[] atomNamesTranslationsMESHI = {
            "AC", "ACA", "ACB", "AN", "AO",
            "CC", "CCA", "CCB", "CN", "CO", "CSG",
            "DC", "DCA", "DCB", "DCG", "DN", "DO", "DOD",
            "EC", "ECA", "ECB", "ECD", "ECG", "EN", "EO", "EOE",
            "FC", "FCA", "FCB", "FCD", "FCE", "FCG", "FCZ", "FN", "FO",
            "GC", "GCA", "GN", "GO",
            "HC", "HCA", "HCB", "HCD", "HCE", "HCG", "HN", "HND", "HNE", "HO",
            "IC", "ICA", "ICB", "ICD", "ICG1", "ICG2", "IN", "IO",
            "KC", "KCA", "KCB", "KCD", "KCE", "KCG", "KN", "KNZ", "KO",
            "LC", "LCA", "LCB", "LCD1", "LCD2", "LCG", "LN", "LO",
            "MC", "MCA", "MCB", "MCE", "MCG", "MN", "MO", "MSD",
            "NC", "NCA", "NCB", "NCG", "NN", "NND", "NO", "NOD",
            "PC", "PCA", "PCB", "PCD", "PCG", "PN", "PO",
            "QC", "QCA", "QCB", "QCD", "QCG", "QN", "QNE", "QO", "QOE",
            "RC", "RCA", "RCB", "RCD", "RCG", "RCZ", "RN", "RNE", "RNH", "RO",
            "SC", "SCA", "SCB", "SN", "SO", "SOG",
            "TC", "TCA", "TCB", "TCG", "TN", "TO", "TOG",
            "VC", "VCA", "VCB", "VCG1", "VCG2", "VN", "VO",
            "WC", "WCA", "WCB", "WCD1", "WCD2", "WCE2", "WCE3", "WCG", "WCH2", "WCZ2", "WCZ3", "WN", "WNE", "WO",
            "YC", "YCA", "YCB", "YCD", "YCE", "YCG", "YCZ", "YN", "YO", "YOH"
    };

    private static String[] atomNamesTranslationsSumma = {
            "AC", "ACA", "ACB", "AN", "AO",
            "CC", "CCA", "CCB", "CN", "CO", "CSG",
            "DC", "DCA", "DCB", "DCG", "DN", "DO", "DOD1",
            "EC", "ECA", "ECB", "ECD", "ECG", "EN", "EO", "EOE1",
            "FC", "FCA", "FCB", "FCD1", "FCE1", "FCG", "FCZ", "FN", "FO",
            "GC", "GCA", "GN", "GO",
            "HC", "HCA", "HCB", "HCD2", "HCE1", "HCG", "HN", "HND1", "HNE2", "HO",
            "IC", "ICA", "ICB", "ICD1", "ICG1", "ICG2", "IN", "IO",
            "KC", "KCA", "KCB", "KCD", "KCE", "KCG", "KN", "KNZ", "KO",
            "LC", "LCA", "LCB", "LCD1", "LCD2", "LCG", "LN", "LO",
            "MC", "MCA", "MCB", "MCE", "MCG", "MN", "MO", "MSD",
            "NC", "NCA", "NCB", "NCG", "NN", "NND2", "NO", "NOD1",
            "PC", "PCA", "PCB", "PCD", "PCG", "PN", "PO",
            "QC", "QCA", "QCB", "QCD", "QCG", "QN", "QNE2", "QO", "QOE1",
            "RC", "RCA", "RCB", "RCD", "RCG", "RCZ", "RN", "RNE", "RNH1", "RO",
            "SC", "SCA", "SCB", "SN", "SO", "SOG",
            "TC", "TCA", "TCB", "TCG2", "TN", "TO", "TOG1",
            "VC", "VCA", "VCB", "VCG1", "VCG2", "VN", "VO",
            "WC", "WCA", "WCB", "WCD1", "WCD2", "WCE2", "WCE3", "WCG", "WCH2", "WCZ2", "WCZ3", "WN", "WNE1", "WO",
            "YC", "YCA", "YCB", "YCD1", "YCE1", "YCG", "YCZ", "YN", "YO", "YOH"
    };

}
