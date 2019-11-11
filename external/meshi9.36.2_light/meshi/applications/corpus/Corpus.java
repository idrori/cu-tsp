/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.*;
import meshi.energy.solvateRot1.SolvateRot1Creator;
import meshi.geometry.DistanceMatrix;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.MeshiPotential;
import meshi.util.Command;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.UpdateableException;
import meshi.util.overlap.Overlap;
import meshi.util.rotamericTools.RotamericTools;

import java.io.*;
import java.text.DecimalFormat;
import java.util.StringTokenizer;

/**
 * A class that hold the 20-type-mutation energy data of one or several proteins.
 * This class can be loaded from disk or created from a PDB file of a protein.
 * Another way to create this class is by merging two instances into a larger corpus.
 */

public class Corpus implements CompositeTorsionsDefinitions, MeshiPotential, KeyWords {

    protected String[] proteinNames = new String[0];
    protected String[] energyNames = new String[0];
    protected int Nenergies = 0;
    protected int Nres = 0;
    protected int[] protInd = new int[0];
    protected int[] resNum = new int[0];
    protected int[] resType = new int[0];
    protected double[] prePro = new double[0];
    protected double[][][] energies = new double[0][][]; // Indexing-[residue index][energy index][mutatation type]
    protected double[][][] coors = new double[0][][]; // Indexing-[residue index][atom: N=0 ; Ca=1 ; C=2][coordinates x=0;y=1;z=2]
    protected double[][] torsions = new double[0][]; // Indexing-[residue index][torsion: omg=0, phi=1, psi=2]
    protected boolean[] resHasAllEne = new boolean[0];
    protected int[] ungapped = null;
    protected String[] excludedProteinsFromUngapped = null;

    protected double[][] forAlign1 = null;
    protected double[][] forAlign2 = null;


    public Corpus(String PDBfile, CommandList commands, EnergyCreator[] energyCreators) throws UpdateableException, EvaluationException{
        proteinNames = new String[1];
        proteinNames[0] = PDBfile;
        Nenergies = energyCreators.length + 2;
        energyNames = new String[Nenergies];
        energyNames[0] = "prop2D";
        energyNames[1] = "Rot1 Solvate";


        // Protein and distance matrix
        Protein protein = new Protein(new AtomList(PDBfile), ResidueExtendedAtomsCreator.creator);
        protein.defrost();

        // finding how many residues with full backbone this protein has
        for (int c = protein.chain().firstNonDummyResidue().number(); c < protein.residues().size(); c++)
            if (haveAllBackboneAtoms(protein.residues().get(c)))
                Nres++;

        // Initializing arrays
        protInd = new int[Nres]; // All the residues come from protein index 0
        resHasAllEne = new boolean[Nres]; // at first the residues are incomplete
        for (int c = 0; c < Nres; c++)
            resHasAllEne[c] = true;
        resNum = new int[Nres];
        resType = new int[Nres];
        energies = new double[Nres][Nenergies][20];
        coors = new double[Nres][3][3];
        torsions = new double[Nres][3];
        prePro = new double[Nres];

        // Propensity parameters
        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        String cpplFileName = parametersDirectory + "/" + COMPOSITE_PROPENSITY_2D_WITH_PP_PARAMETERS;
        SplinedPolynomialsLoader spl;
        try {
            spl = new SplinedPolynomialsLoader(cpplFileName);
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
        SplinedPolynomial[] splines = new SplinedPolynomial[22];
        for (int mutateTo = 0; mutateTo < 22; mutateTo++)
            splines[mutateTo] = spl.findPolynomial(mutateTo, POLYNOMIAL_PHI_PSI_TORSIONS, ALL);

        // Running isOn the residues extracting structural data and propensity
        Nres = 0;
        ResidueTorsionsList rtl = new ResidueTorsionsList(protein, protein.atoms().molecularSystem().getDistanceMatrix());
        for (int c = protein.chain().firstNonDummyResidue().number(); c < protein.residues().size(); c++)
            if (haveAllBackboneAtoms(protein.residues().get(c))) {
                Residue res = protein.residues().get(c);
                resType[Nres] = protein.residues().get(c).type.ordinal();
                resNum[Nres] = protein.residues().get(c).number();
                extractCoordinatesOfRes(res, coors[Nres]);
                ResidueTorsions rt = rtl.findResidueInList(res.number());
                if (rt != null) {
                    if (rt.getTorsion(OMG) != null)
                        torsions[Nres][0] = rt.getTorsion(OMG).torsion();
                    else
                        torsions[Nres][0] = Double.NaN;
                    if (rt.getTorsion(PHI) != null)
                        torsions[Nres][1] = rt.getTorsion(PHI).torsion();
                    else
                        torsions[Nres][1] = Double.NaN;
                    if (rt.getTorsion(PSI) != null)
                        torsions[Nres][2] = rt.getTorsion(PSI).torsion();
                    else
                        torsions[Nres][2] = Double.NaN;
                }
                // Calculating the prop2D energies
                if ((rt.getTorsion(PSI) != null) && (rt.getTorsion(PHI) != null)) {
                    for (int mutateTo = 0; mutateTo < 20; mutateTo++)
                        energies[Nres][0][mutateTo] = splines[mutateTo].value(0, torsions[Nres][1], torsions[Nres][2]) -
                                splines[OMNI].value(0, torsions[Nres][1], torsions[Nres][2]) + 1e-100;
                    prePro[Nres] = splines[PREPRO].value(0, torsions[Nres][1], torsions[Nres][2]) -
                            splines[OMNI].value(0, torsions[Nres][1], torsions[Nres][2]) + 1e-100;
                }
                Nres++;
            }

        // Calculating Solvate Rot1
        // ************************
        calculateSolRot1(PDBfile, commands);

        // Checking if all residues have Kosher energies
        for (int Eind = 0; Eind < Nenergies; Eind++) {
            for (int c = 0; c < Nres; c++)
                for (int mutateTo = 0; mutateTo < 20; mutateTo++)
                    if (energies[c][Eind][mutateTo] == 0.0)
                        resHasAllEne[c] = false;
        }
    }

    /**
     * Reading a corpus from file. See the file format in the method 'writeToDisk'
     */
    public Corpus(String exsitingCorpusFile) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(exsitingCorpusFile));
            String line = br.readLine();
            Nenergies = (new Integer(line)).intValue();
            energyNames = new String[Nenergies];
            for (int Eind = 0; Eind < Nenergies; Eind++)
                energyNames[Eind] = br.readLine().trim();
            line = br.readLine();
            proteinNames = new String[(new Integer(line)).intValue()];
            for (int nProt = 0; nProt < proteinNames.length; nProt++)
                proteinNames[nProt] = br.readLine().trim();
            line = br.readLine();
            Nres = (new Integer(line)).intValue();
            protInd = new int[Nres];
            resHasAllEne = new boolean[Nres];
            for (int c = 0; c < Nres; c++)
                resHasAllEne[c] = true;
            resNum = new int[Nres];
            resType = new int[Nres];
            energies = new double[Nres][Nenergies][20];
            coors = new double[Nres][3][3];
            torsions = new double[Nres][3];
            prePro = new double[Nres];
            for (int c = 0; c < Nres; c++) {
                StringTokenizer st = new StringTokenizer(br.readLine());
                protInd[c] = (new Integer(st.nextToken())).intValue();
                resNum[c] = (new Integer(st.nextToken())).intValue();
                resType[c] = (new Integer(st.nextToken())).intValue();
                for (int c1 = 0; c1 < 3; c1++)
                    for (int c2 = 0; c2 < 3; c2++)
                        coors[c][c1][c2] = (new Double(st.nextToken())).doubleValue();
                for (int c1 = 0; c1 < 3; c1++)
                    torsions[c][c1] = (new Double(st.nextToken())).doubleValue();
                prePro[c] = (new Double(st.nextToken())).doubleValue();
                for (int Eind = 0; Eind < Nenergies; Eind++) {
                    st = new StringTokenizer(br.readLine());
                    for (int mutateTo = 0; mutateTo < 20; mutateTo++)
                        energies[c][Eind][mutateTo] = (new Double(st.nextToken())).doubleValue();
                }
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Writing all the residues in the corpus that are marked as TRUE in the 'toWrite' input array (its length must be as the
     * corpus's). The format is:
     * <p/>
     * <num of energies>
     * <Energy 1 - name>
     * <Energy 2 - name>
     * ...
     * <Nenergy Nenergy - name>
     * <Num of Proteins>
     * <Protein 1 - name>
     * <Protein 2 - name>
     * ...
     * <Protein Nprot - name>
     * <Num of residues>
     * <residue 1 - <prot index> <res num> <res type> <Nx> <Ny> <Nz> <CAx> <CAy> <CAz> <Cx> <Cy> <Cz> <OMG> <PHI> <PSI> <PrePro>>
     * <residue 1 - <Energy 1 in ALA> <Energy 1 in CYS> ... <Energy 1 in TYR> <Energy 2 in ALA> <Energy 2 in CYS> ... <Energy 2 in TYR> ... >
     * <residue 2 - <prot index> <res num> <res type> <Nx> <Ny> <Nz> <CAx> <CAy> <CAz> <Cx> <Cy> <Cz> <OMG> <PHI> <PSI> <PrePro>>
     * <residue 2 - <Energy 1 in ALA> <Energy 1 in CYS> ... <Energy 1 in TYR> <Energy 2 in ALA> <Energy 2 in CYS> ... <Energy 2 in TYR> ... >
     * ....
     * <residue Nres - <prot index> <res num> <res type> <Nx> <Ny> <Nz> <CAx> <CAy> <CAz> <Cx> <Cy> <Cz> <OMG> <PHI> <PSI> <PrePro>>
     * <residue Nres - <Energy 1 in ALA> <Energy 1 in CYS> ... <Energy 1 in TYR> <Energy 2 in ALA> <Energy 2 in CYS> ... <Energy 2 in TYR> ... >
     */
    private void writeToDisk(String fileName, boolean[] toWrite) {
        DecimalFormat fmt = new DecimalFormat("0.####");
        int TrueNres = 0;
        for (int c = 0; c < Nres; c++)
            if (toWrite[c])
                TrueNres++;
        try {
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(fileName)));
            pw.println(Nenergies);
            for (int Eind = 0; Eind < Nenergies; Eind++)
                pw.println(energyNames[Eind]);
            pw.println(proteinNames.length);
            for (int nProt = 0; nProt < proteinNames.length; nProt++)
                pw.println(proteinNames[nProt]);
            pw.println(TrueNres);
            for (int c = 0; c < Nres; c++)
                if (toWrite[c]) {
                    pw.print(protInd[c] + " " + resNum[c] + " " + resType[c] + " ");
                    for (int c1 = 0; c1 < 3; c1++)
                        for (int c2 = 0; c2 < 3; c2++)
                            pw.print(fmt.format(coors[c][c1][c2]) + " ");
                    for (int c1 = 0; c1 < 3; c1++)
                        pw.print(fmt.format(torsions[c][c1]) + " ");
                    pw.print(prePro[c]);
                    pw.println();
                    for (int Eind = 0; Eind < Nenergies; Eind++) {
                        for (int mutateTo = 0; mutateTo < 20; mutateTo++)
                            pw.print(energies[c][Eind][mutateTo] + " ");
                        pw.println();
                    }
                }
            pw.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * This overload writes the entire corpus to disk. Residues that do not have all the energies are not written.
     */
    public void writeToDisk(String fileName) {
        boolean[] toWrite = new boolean[Nres];
        for (int c = 0; c < Nres; c++)
            if (resHasAllEne[c])
                toWrite[c] = true;
            else
                toWrite[c] = false;
        writeToDisk(fileName, toWrite);
    }

    /**
     * This overload writes to disk a subset of the corpus, all the residues that have of the previous method, that writes to disk only residues from
     * the fragments whos indices are given in the input array. This is useful to create
     * more compact libraries from a big one. The format is as for the other method.
     */
    public void writeToDisk(String fileName, int indAr, int fragL) {
    }

    public void merge(Corpus corpus) {
        // Checking compatibility
        if (Nenergies != corpus.Nenergies)
            throw new RuntimeException("Can not merge corpura. Number of energies mismatch");
        for (int Eind = 0; Eind < Nenergies; Eind++)
            if (!energyNames[Eind].equals(corpus.energyNames[Eind]))
                throw new RuntimeException("Can not merge corpura. Energy mismatch at energy No. " + (Eind + 1));

        // Merging
        ungapped = null;
        String[] tmp_proteinNames = new String[proteinNames.length + corpus.proteinNames.length];
        int tmp_Nres = Nres + corpus.Nres;
        int[] tmp_protInd = new int[protInd.length + corpus.protInd.length];
        int[] tmp_resNum = new int[resNum.length + corpus.resNum.length];
        int[] tmp_resType = new int[resType.length + corpus.resType.length];
        double[][][] tmp_energies = new double[energies.length + corpus.energies.length][][];
        double[][][] tmp_coors = new double[coors.length + corpus.coors.length][][];
        double[][] tmp_torsions = new double[torsions.length + corpus.torsions.length][];
        double[] tmp_prePro = new double[prePro.length + corpus.prePro.length];
        boolean[] tmp_resHasAllEne = new boolean[resHasAllEne.length + corpus.resHasAllEne.length];

        for (int c = 0; c < proteinNames.length; c++)
            tmp_proteinNames[c] = proteinNames[c];
        for (int c = 0; c < corpus.proteinNames.length; c++)
            tmp_proteinNames[c + proteinNames.length] = corpus.proteinNames[c];

        for (int c = 0; c < protInd.length; c++)
            tmp_protInd[c] = protInd[c];
        for (int c = 0; c < corpus.protInd.length; c++)
            tmp_protInd[c + protInd.length] = corpus.protInd[c] + proteinNames.length;

        for (int c = 0; c < resNum.length; c++)
            tmp_resNum[c] = resNum[c];
        for (int c = 0; c < corpus.resNum.length; c++)
            tmp_resNum[c + resNum.length] = corpus.resNum[c];

        for (int c = 0; c < resType.length; c++)
            tmp_resType[c] = resType[c];
        for (int c = 0; c < corpus.resType.length; c++)
            tmp_resType[c + resType.length] = corpus.resType[c];

        for (int c = 0; c < prePro.length; c++)
            tmp_prePro[c] = prePro[c];
        for (int c = 0; c < corpus.prePro.length; c++)
            tmp_prePro[c + prePro.length] = corpus.prePro[c];

        for (int c = 0; c < energies.length; c++)
            tmp_energies[c] = energies[c];
        for (int c = 0; c < corpus.energies.length; c++)
            tmp_energies[c + energies.length] = corpus.energies[c];

        for (int c = 0; c < coors.length; c++)
            tmp_coors[c] = coors[c];
        for (int c = 0; c < corpus.coors.length; c++)
            tmp_coors[c + coors.length] = corpus.coors[c];

        for (int c = 0; c < torsions.length; c++)
            tmp_torsions[c] = torsions[c];
        for (int c = 0; c < corpus.torsions.length; c++)
            tmp_torsions[c + torsions.length] = corpus.torsions[c];

        for (int c = 0; c < resHasAllEne.length; c++)
            tmp_resHasAllEne[c] = resHasAllEne[c];
        for (int c = 0; c < corpus.resHasAllEne.length; c++)
            tmp_resHasAllEne[c + resHasAllEne.length] = corpus.resHasAllEne[c];

        proteinNames = tmp_proteinNames;
        Nres = tmp_Nres;
        protInd = tmp_protInd;
        resNum = tmp_resNum;
        resType = tmp_resType;
        energies = tmp_energies;
        coors = tmp_coors;
        torsions = tmp_torsions;
        resHasAllEne = tmp_resHasAllEne;
        prePro = tmp_prePro;
    }


    public void buildUngappedArray(int fragL) {
        buildUngappedArray(fragL, excludedProteinsFromUngapped);
    }

    /**
     * Building the ungapped array for fragments of length fragL.
     * Fragments from the proteins in the provided list are excluded.
     */
    public void buildUngappedArray(int fragL, String[] excludeProteins) {
        int Nfrag = 0;
        if (excludeProteins == null)
            excludeProteins = new String[0];
        boolean inExcluded = false;
        for (int c = 0; c < Nres - fragL; c++) {
            inExcluded = false;
            for (int d = 0; d < excludeProteins.length; d++)
                if (proteinNames[protInd[c]].indexOf(excludeProteins[d]) > -1)
                    inExcluded = true;
            if (!inExcluded && (protInd[c] == protInd[c + fragL - 1]) && ((resNum[c] + fragL - 1) == resNum[c + fragL - 1]))
                Nfrag++;
        }
        ungapped = new int[Nfrag];
        Nfrag = 0;
        for (int c = 0; c < Nres - fragL; c++) {
            inExcluded = false;
            for (int d = 0; d < excludeProteins.length; d++)
                if (proteinNames[protInd[c]].indexOf(excludeProteins[d]) > -1)
                    inExcluded = true;
            if (!inExcluded && (protInd[c] == protInd[c + fragL - 1]) && ((resNum[c] + fragL - 1) == resNum[c + fragL - 1])) {
                ungapped[Nfrag] = c;
                Nfrag++;
            }
        }

        // Making instances to save time in RMS calculations
        forAlign1 = new double[3][fragL * 3];
        forAlign2 = new double[3][fragL * 3];
    }

    public void threadingExperiment(int fragL) {
        double nat = 0, thr = 0;
        int[] seq = new int[fragL];
        int rank = 0;
        buildUngappedArray(fragL);
        System.out.println("***********************");
        System.out.println("Number of fragments = " + ungapped.length);
        System.out.println("***********************");
        for (int c = 0; c < ungapped.length; c += 10) {
            for (int e = 0; e < fragL; e++)
                seq[e] = resType[ungapped[c] + e];
            nat = 0;
            for (int e = 0; e < fragL; e++)
                nat += energies[ungapped[c] + e][0][seq[e]];
            rank = 0;
            for (int d = 0; d < ungapped.length; d++) {
                thr = 0;
                for (int e = 0; e < fragL; e++)
                    thr += energies[ungapped[d] + e][0][seq[e]];
                if (thr <= nat)
                    rank++;
            }
            System.out.println(rank);
        }
    }

    public void threadingExperiment_withPP(int fragL, double Wprop, double Wsolv) {
        double nat = 0, thr = 0;
        int[] seq = new int[fragL];
        boolean[] pp = new boolean[fragL];
        for (int e = 0; e < fragL; e++)
            pp[e] = false;
        int rank = 0;
        buildUngappedArray(fragL);
        System.out.println("***********************");
        System.out.println("Number of fragments = " + ungapped.length);
        System.out.println("***********************");
        for (int c = 0; c < ungapped.length; c += 10) {
            for (int e = 0; e < fragL; e++)
                seq[e] = resType[ungapped[c] + e];
            for (int e = 0; e < fragL - 1; e++)
                if ((seq[e + 1] == 12) && (seq[e] != 5) && (seq[e] != 12))
                    pp[e] = true;
                else
                    pp[e] = false;
            pp[fragL - 1] = false;
            nat = 0;
            for (int e = 0; e < fragL; e++) {
                if (pp[e])
                    nat += Wprop * prePro[ungapped[c] + e];
                else
                    nat += Wprop * energies[ungapped[c] + e][0][seq[e]];
                nat += Wsolv * energies[ungapped[c] + e][1][seq[e]];
            }
            rank = 0;
            for (int d = 0; d < ungapped.length; d++) {
                thr = 0;
                for (int e = 0; e < fragL; e++) {
                    if (pp[e])
                        thr += Wprop * prePro[ungapped[d] + e];
                    else
                        thr += Wprop * energies[ungapped[d] + e][0][seq[e]];
                    thr += Wsolv * energies[ungapped[d] + e][1][seq[e]];
                }
                if (thr <= nat)
                    rank++;
            }
            System.out.println(rank);
        }
    }

    /**
     * This is a method for general ungapped threading. You give it the frag length 'fragL', and the number of instances
     * you want 'Ninstances'. It will randomly pick 'Ninstances' random fragments and thread the sequence of the first into
     * all the others. Than print:
     * <header> <1> <RMS to frag 1> <frag 1 propensity> <frag 1 blosum62>
     * <header> <2> <RMS to frag 1> <frag 2 propensity> <frag 2 blosum62>
     * <header> <3> <RMS to frag 1> <frag 3 propensity> <frag 3 blosum62>
     * ...
     * <header> <Ninstances> <RMS to frag 1> <frag Ninstances propensity> <frag Ninstances blosum62>
     * <p/>
     * The RMS is calculated between two indices (inclusive) in the FRAG reference frame.
     */

    public void GeneralThreadingExperiment(int fragL, int Ninstances, String header, int RMSind1, int RMSind2) {
        if (ungapped == null)
            buildUngappedArray(fragL);
        // the random first fragment
        int natInd = (int) (ungapped.length * Math.random());
        // Finding the sequence and pre-PRO
        int[] seq = new int[fragL];
        boolean[] pp = new boolean[fragL];
        for (int e = 0; e < fragL; e++)
            pp[e] = false;
        for (int e = 0; e < fragL; e++)
            seq[e] = resType[ungapped[natInd] + e];
        for (int e = 0; e < fragL - 1; e++)
            if ((seq[e + 1] == 12) && (seq[e] != 5) && (seq[e] != 12))
                pp[e] = true;
            else
                pp[e] = false;
        pp[fragL - 1] = false;
        // Looping isOn the instances
        for (int cc = 0; cc < Ninstances; cc++) {
            int randomPick;
            if (cc == 0)
                randomPick = natInd;
            else
                randomPick = (int) (ungapped.length * Math.random());
            // Printing the RMS
            System.out.print(header + " " + (cc + 1) + " " +
                    calcRmsBetweenStruct(ungapped[natInd] + RMSind1, ungapped[randomPick] + RMSind1, RMSind2 - RMSind1 + 1, 1, RMSind2 - RMSind1 + 1)
                    + " ");
            // Threading
            double propScore = 0.0;
            double blosumScore = 0.0;
            for (int e = 0; e < fragL; e++) {
                if (pp[e])
                    propScore += prePro[ungapped[randomPick] + e];
                else
                    propScore += energies[ungapped[randomPick] + e][0][seq[e]];
                blosumScore -= blosum62[seq[e]][resType[ungapped[randomPick] + e]];
            }
            System.out.println(propScore + " " + blosumScore);
        }
    }

// Auxilary method

    private static boolean haveAllBackboneAtoms(Residue res) {
        return ((res.ca() != null) &&
                (!res.ca().nowhere()) &&
                (res.amideN() != null) &&
                (!res.amideN().nowhere()) &&
                (res.carbonylC() != null) &&
                (!res.carbonylC().nowhere()));
    }

// Auxilary method

    private static void extractCoordinatesOfRes(Residue res, double[][] co) {
        Atom ca = res.ca();
        Atom amideN = res.amideN();
        Atom carbonylC = res.carbonylC();
        co[0][0] = amideN.x();    //res.getAtoms().findAtomInList("N" , res.number).x();
        co[0][1] = amideN.y();    //res.getAtoms().findAtomInList("N" , res.number).y();
        co[0][2] = amideN.z();    //res.getAtoms().findAtomInList("N" , res.number).z();
        co[1][0] = ca.x();        //res.getAtoms().findAtomInList("CA" , res.number).x();
        co[1][1] = ca.y();        //res.getAtoms().findAtomInList("CA" , res.number).y();
        co[1][2] = ca.z();        //res.getAtoms().findAtomInList("CA" , res.number).z();
        co[2][0] = carbonylC.x(); //res.getAtoms().findAtomInList("C" , res.number).x();
        co[2][1] = carbonylC.y(); //res.getAtoms().findAtomInList("C" , res.number).y();
        co[2][2] = carbonylC.z(); //res.getAtoms().findAtomInList("C" , res.number).z();
    }


// Method for calculating Solvate Rot1 in the 20 mutations
// ************************

    private void calculateSolRot1(String PDBfile, CommandList commands) throws UpdateableException,EvaluationException{
        // making the auxilary protein, and putting it in Rot1
        AtomList al1 = new AtomList(PDBfile);
        al1.moveCMtoOrigin();
        AtomList al2 = new AtomList("tmpPdb.pdb");
        for (int c = 0; c < al2.size(); c++)
            al1.add(al2.atomAt(c));
        Protein aux = new Protein(al1, ResidueExtendedAtomsCreator.creator);
        aux.defrost();
        DistanceMatrix dm = aux.atoms().molecularSystem().getDistanceMatrix();

        double[][] allPP = RotamericTools.phipsi(aux, dm);
        DunbrackLib lib = new DunbrackLib(commands, 1.0, 2);
        double[][] pp = RotamericTools.putIntoRot1(aux, dm, lib);

        for (int c = 1; c < (Nres - 1); c++)
            if ((resNum[c - 1] == resNum[c] - 1) && (resNum[c + 1] == resNum[c] + 1)) { // This residue is not at the termini
                int num = resNum[c];
                System.out.println("AAAAAA: " + num);
                for (int d = 1005; d <= 1059; d += 3)    // No Alanine
                    if (d != 1017) {                       // No glycines
                        pp[d][0] = allPP[num][0];
                        pp[d][1] = allPP[num][1];
                        pp[d][2] = (d - 1002) / 3;
                        pp[d][3] = lib.getRotamerProb((int) pp[d][2], pp[d][0], pp[d][1], 0);
                    }
                EnergyCreator[] energyCreators = {new SolvateRot1Creator(1.0, pp, 1000)};
                dm = aux.atoms().molecularSystem().getDistanceMatrix();
                TotalEnergy energy = new TotalEnergy(aux, energyCreators, commands,"Generic total energy");
                int lastMovedRes = num;
                energies[c][1][0] = energies[c][1][5] = 1e-100;
                for (int mutateTo = 1; mutateTo < 20; mutateTo++)      // No Alanine
                    if (mutateTo != 5) {                                   // No glycines
                        specialMutate(aux, lastMovedRes, 1002 + 3 * mutateTo, allPP, lib);
                        lastMovedRes = 1002 + 3 * mutateTo;
                        try {
                            energy.evaluate();
                            energy.resetAtomEnergies();
                            energy.evaluateAtoms();
                        } catch (UpdateableException ae) {
                            throw new RuntimeException(" calculateSolRot failed due to " + ae);
                        }
                        energies[c][1][mutateTo] = aux.residue(1002 + 3 * mutateTo).ca().energy() + 1e-100;
                    }
                specialMutate(aux, lastMovedRes, num, allPP, lib);
            }
        // aux.getAtoms().print();
    }


    // This auxilary methods:
    // Move the coordinates of mutateTo and the flanks backbone, and put them isOn mutate this.
    // Take the coordinates of mutateThis and the flanks backbone, and add 20 to all their getAtoms.
    // Build the sidechain of mutateTo into Rot1.

    private void specialMutate(Protein aux, int mutateThis, int mutateTo, double[][] pp, DunbrackLib lib) {
        Atom atom;
        atom = aux.residue(mutateThis - 1).amideN();                          //aux.residue(mutateThis-1).getAtoms().findAtomInList("N" , mutateThis-1);
        aux.residue(mutateTo - 1).amideN().setXYZ(atom.x(), atom.y(), atom.z());//aux.residue(mutateTo-1).getAtoms().findAtomInList("N",mutateTo-1).setXYZ(atom.x(),atom.y(),atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis - 1).ca();                          //aux.residue(mutateThis-1).getAtoms().findAtomInList("CA" , mutateThis-1);
        aux.residue(mutateTo - 1).ca().setXYZ(atom.x(), atom.y(), atom.z());//aux.residue(mutateTo-1).getAtoms().findAtomInList("CA" , mutateTo-1).setXYZ(atom.x(),atom.y(),atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis - 1).carbonylC();                           //aux.residue(mutateThis-1).getAtoms().findAtomInList("C" , mutateThis-1);
        aux.residue(mutateTo - 1).carbonylC().setXYZ(atom.x(), atom.y(), atom.z()); //aux.residue(mutateTo-1).getAtoms().findAtomInList("C",mutateTo-1).setXYZ(atom.x(),atom.y(),atom.z(
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis - 1).carbonylO();
        aux.residue(mutateTo - 1).carbonylO().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);


        atom = aux.residue(mutateThis).amideN();
        aux.residue(mutateTo).amideN().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis).ca();
        aux.residue(mutateTo).ca().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis).carbonylC();
        aux.residue(mutateTo).carbonylC().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis).carbonylO();
        aux.residue(mutateTo).carbonylO().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);


        atom = aux.residue(mutateThis + 1).amideN();
        aux.residue(mutateTo + 1).amideN().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis + 1).ca();
        aux.residue(mutateTo + 1).ca().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis + 1).carbonylC();
        aux.residue(mutateTo + 1).carbonylC().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        atom = aux.residue(mutateThis + 1).carbonylO();
        aux.residue(mutateTo + 1).carbonylO().setXYZ(atom.x(), atom.y(), atom.z());
        atom.setXYZ(atom.x() + 20, atom.y() + 20, atom.z() + 20);

        if ((((int) pp[mutateTo][2]) != 0) && (((int) pp[mutateTo][2]) != 5))
            ResidueBuilder.build(aux.residue(mutateTo), (int) pp[mutateTo][2],
                    lib.getRotamer((int) pp[mutateTo][2], pp[mutateTo][0], pp[mutateTo][1], 0));
        else
            ResidueBuilder.build(aux.residue(mutateTo), (int) pp[mutateTo][2], null);
        if ((((int) pp[mutateThis][2]) != 0) && (((int) pp[mutateThis][2]) != 5))
            ResidueBuilder.build(aux.residue(mutateThis), (int) pp[mutateThis][2],
                    lib.getRotamer((int) pp[mutateThis][2], pp[mutateThis][0], pp[mutateThis][1], 0));
        else
            ResidueBuilder.build(aux.residue(mutateThis), (int) pp[mutateThis][2], null);
    }


    public int findInUngapped(String proteinName, int residueNumber) {
        proteinName = proteinName.toUpperCase().trim();
        for (int c = 0; c < ungapped.length; c++)
            if ((proteinNames[protInd[ungapped[c]]].toUpperCase().indexOf(proteinName) > -1) &&
                    (resNum[ungapped[c]] == residueNumber))
                return c;
        return -1;
    }

    /**
     * This method gives the rms between two fragments in the corpus (starting in ind1 and ind2). fragL is the
     * fragment length. The fragments are superimposed according to an overlap region who's size is a parameter. The type of
     * the overlap region are given by manner:
     * -1  - overlap of the N-terminus
     * 0 - two overlap regions in both termini (each of length overlap)
     * 1  - overlap of the C-terminus
     * You can force the overlap of the entire fragments by giving (overlap=fragL) and (manner=1)
     */
    public double calcRmsBetweenStruct(int ind1, int ind2, int fragL, int manner, int overlap) {
        if (forAlign1[0].length != fragL) {
            forAlign1 = new double[3][fragL * 3];
            forAlign2 = new double[3][fragL * 3];
        }
        // Making dummy 1 ready
        for (int c = 0; c < fragL; c++)
            for (int c1 = 0; c1 < 3; c1++)
                for (int c2 = 0; c2 < 3; c2++)
                    forAlign1[c2][c * 3 + c1] = coors[ind1 + c][c1][c2];
        // Making dummy 2 ready
        for (int c = 0; c < fragL; c++)
            for (int c1 = 0; c1 < 3; c1++)
                for (int c2 = 0; c2 < 3; c2++)
                    forAlign2[c2][c * 3 + c1] = coors[ind2 + c][c1][c2];
        //Overlap overlap = new Overlap(dummy1, dummy2, fragL*3, "1", "2");
        //return overlap.rms();
        if (manner == 0) {   // Both sides
            int[] ar = new int[2 * 3 * overlap];
            for (int c = 0; c < 3 * overlap; c++) {
                ar[c] = c;
                ar[2 * 3 * overlap - 1 - c] = 3 * fragL - 1 - c;
            }
            try {
                return Overlap.rmsPartial(forAlign1, forAlign2, ar);
            }
            catch (Exception e) {
                return -1;
            }
        } else if (manner == 1) {   // C terminus
            int[] ar = new int[3 * overlap];
            for (int c = 0; c < 3 * overlap; c++) {
                ar[3 * overlap - 1 - c] = 3 * fragL - 1 - c;
            }
            try {
                return Overlap.rmsPartial(forAlign1, forAlign2, ar);
            }
            catch (Exception e) {
                return -1;
            }
        } else if (manner == -1) {   // N terminus
            int[] ar = new int[3 * overlap];
            for (int c = 0; c < 3 * overlap; c++) {
                ar[c] = c;
            }
            try {
                return Overlap.rmsPartial(forAlign1, forAlign2, ar);
            }
            catch (Exception e) {
                return -1;
            }
        } else
            throw new RuntimeException("manner parameter must be -1,0 or +1");
    }


    public void buildUngappedArraynoPG(int fragL) {
        int Nfrag = 0;
        for (int c = 0; c < Nres - fragL; c++)
            if ((protInd[c] == protInd[c + fragL - 1]) && ((resNum[c] + fragL - 1) == resNum[c + fragL - 1])) {
                boolean withpg = false;
                for (int d = 0; d < fragL; d++)
                    if ((resType[c + d] == 5) || (resType[c + d] == 12))
                        withpg = true;
                if (!withpg)
                    Nfrag++;
            }
        ungapped = new int[Nfrag];
        Nfrag = 0;
        for (int c = 0; c < Nres - fragL; c++)
            if ((protInd[c] == protInd[c + fragL - 1]) && ((resNum[c] + fragL - 1) == resNum[c + fragL - 1])) {
                boolean withpg = false;
                for (int d = 0; d < fragL; d++)
                    if ((resType[c + d] == 5) || (resType[c + d] == 12))
                        withpg = true;
                if (!withpg) {
                    ungapped[Nfrag] = c;
                    Nfrag++;
                }
            }

        // Making instances to save time in RMS calculations
        forAlign1 = new double[3][fragL * 3];
        forAlign2 = new double[3][fragL * 3];
    }

    public void threadingExperimentnoPG(int fragL) {
        double nat = 0, thr = 0;
        int[] seq = new int[fragL];
        int rank = 0;
        buildUngappedArraynoPG(fragL);
        System.out.println("***********************");
        System.out.println("Number of fragments = " + ungapped.length);
        System.out.println("***********************");
        for (int c = 0; c < ungapped.length; c += 10) {
            for (int e = 0; e < fragL; e++)
                seq[e] = resType[ungapped[c] + e];
            nat = 0;
            for (int e = 0; e < fragL; e++)
                nat += energies[ungapped[c] + e][0][seq[e]];
            rank = 0;
            for (int d = 0; d < ungapped.length; d++) {
                thr = 0;
                for (int e = 0; e < fragL; e++)
                    thr += energies[ungapped[d] + e][0][seq[e]];
                if (thr <= nat)
                    rank++;
            }
            System.out.println(rank);
        }
    }


    /**
     * Find the distributaion of distances between ends.
     */
    public void findDisDist(int len) {
        DecimalFormat fmt = new DecimalFormat("0.#");
        buildUngappedArray(len);
        //dists = new double[ungapped.length];
        for (int c = 0; c < 50000; c++) { //ungapped.length ; c++) {
            System.out.print(fmt.format(Math.sqrt((coors[ungapped[c]][0][0] - coors[ungapped[c] + len - 1][2][0]) *
                    (coors[ungapped[c]][0][0] - coors[ungapped[c] + len - 1][2][0]) +
                    (coors[ungapped[c]][0][1] - coors[ungapped[c] + len - 1][2][1]) *
                            (coors[ungapped[c]][0][1] - coors[ungapped[c] + len - 1][2][1]) +
                    (coors[ungapped[c]][0][2] - coors[ungapped[c] + len - 1][2][2]) *
                            (coors[ungapped[c]][0][2] - coors[ungapped[c] + len - 1][2][2]))) + " ");
        }
    }

    /**
     * Find the relation between Delta{phi,psi} and RMS
     */
    public void findPhiPsiRelations(int howMany) {
        buildUngappedArray(2);
        for (int c = 0; c < howMany; c++) {
            int ind1 = (int) (ungapped.length * Math.random());
            int ind2 = (int) (ungapped.length * Math.random());
            double delta1 = Math.abs(torsions[ungapped[ind1]][2] - torsions[ungapped[ind2]][2]);
            double delta2 = Math.abs(torsions[ungapped[ind1] + 1][1] - torsions[ungapped[ind2] + 1][1]);
            if (delta1 > Math.PI)
                delta1 = 2 * Math.PI - delta1;
            if (delta2 > Math.PI)
                delta2 = 2 * Math.PI - delta2;
            System.out.print((delta1 * delta1 + delta2 * delta2) + " ");
            System.out.print(calcRmsBetweenStruct(ungapped[ind1], ungapped[ind2], 2, 1, 2) + "\n");
        }
    }


    /* BLOSUM62. Indices are A,C,D,E,F,...,W,Y,X */
    public final static double[][] blosum62 =
            {{4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2, 0},
                    {0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2, -2},
                    {-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3, -1},
                    {-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2, -1},
                    {-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3, -1},
                    {0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3, -1},
                    {-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2, -1},
                    {-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1, -1},
                    {-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2, -1},
                    {-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1, -1},
                    {-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1, -1},
                    {-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2, -1},
                    {-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3, -2},
                    {-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1, -1},
                    {-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2, -1},
                    {1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2, 0},
                    {0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2, 0},
                    {0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1, -1},
                    {-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2, -2},
                    {-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7, -1},
                    {0, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1, 0, 0, -1, -2, -1, -1}};

    public void setExcludedProteinsFromUngapped(String[] list) {
        excludedProteinsFromUngapped = list;
    }

} // Of Corpus
