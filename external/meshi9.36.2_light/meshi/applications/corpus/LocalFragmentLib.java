/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import meshi.energy.simpleEnergyTerms.compositeTorsions.CompositeTorsionsDefinitions;
import meshi.geometry.ResidueBuilder;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.util.*;
import meshi.util.overlap.*;

/**
 * The purpose of this class is is class is to build a fragment library that is compatible with a certain
 * sequence. The lib is applicable to a corpus instance. The sequence is given to the constructor as an array
 * of integers (in the range: 0-19). The constructor creates a library with at most libSize fragments that are
 * different from each other by at least a certain RMS curoff. The final size of the lib could be smaller than
 * the requested if this cutoff is stringent. In general the RMS (for fragment diversity purposes) is calculated
 * not isOn the entire sequence, but rather isOn an inner subsequence. I found out that when you do threading
 * isOn a long sequence the results are more accurate for the inner subsequences. The ranking of fragments in the lib
 * is currently by propensity energy and pair-wise sequence alignment using BLOSUM62.
 */

public class LocalFragmentLib implements CompositeTorsionsDefinitions, MeshiPotential,
        KeyWords {

    private Corpus corpus = null;
    private int fragL = 0;
    private int libSize = 0;
    private double rmsCutOff;
    private int[] seq = null;
    private int[] libOrig = null;  // This is the lib core. The indices here point to the fragments in the corpus. The size of this array is the lib size.
    private double[] libEnergy = null;  // This is the energy ascociated with a fragment. The size of this array is the lib size.
    private boolean[] similarToLib = null; // This marks (true values) fragments that are already selected, or are very similar to fragments alreay in the lib. The size of this array is as curpus.ungapped
    private int trueFragStart = -999;  // The place in the truncated (no extension) sequence where the true frag (no overlap) starts
    private int trueFragEnd = -999;  // The place in the truncated (no extension) sequence where the true frag (no overlap) ends


//	These fields are needed when the library is constructed so that it is compatible with the stalks.
    private Protein prot = null;
    private int residueFragStart = -999;
    private int overlap = -999;
    private int manner = -999;
    private double overlapRMSsimilarity = 0.0;

    private double[][] fragCoors = null; // This is an array that is used a lot in: 'isFragCompatible'. It is the coordinates of the fragment
    private double[][] protCoors = null; // This is an array that is used a lot in: 'isFragCompatible'. It is the coordinates of the protein


    /* This constructor is able to construct a frag-lib, with some lateral addendums isOn the fragment as follows:
      *
      * |-------------------------------------------SEQ------------------------------------------------------------------------|
      * [   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ][   ]
      *                          |--overlap----||--------------------frag--------------||--overlap----|
      * |---------expandThreadingN-------------|                                       |-----------------expandThreadingC------|
      *
      * Where "frag" is what we are interested in. We make sure, all the "frag" parts are at least 'rmsCutoff' mutually different than
      * all the other "frags" in the library. In this constructor this difference is calculated after the OVERLAP regions
      * are superimposed. Overlap may be isOn one side of "frag" or both, depending isOn 'manner' (which can be -1,0,1). The overlap
      * must also have an RMS of less than 'overlapRMSsimilarity' to the protein in the overlapping region.
      *
      * The parameter 'fragStartInSeq' is the index in 'seq' (begining of seq is 0) where 'frag' or 'overlap' (the first of the two) begin.
      * The parameter 'fragEndInSeq' is the index in 'seq' (begining of seq is 0) where 'frag' or 'overlap' (the latter of the two) end.
      * The parameter 'residueFragStart' should be the residue number in 'prot' where where 'frag' or 'overlap' (the first of the two) begin.
      *
      */

    public LocalFragmentLib(Corpus corpus, int[] seq, int libSize, double rmsCutoff, int fragStartInSeq, int fragEndInSeq,
                            Protein prot, int residueFragStart, int overlap, int manner, double overlapRMSsimilarity) {
        this.corpus = corpus;
        this.seq = seq;
        this.libSize = libSize;
        this.rmsCutOff = rmsCutOff;
        this.prot = prot;
        this.residueFragStart = residueFragStart;
        this.overlap = overlap;
        this.manner = manner;
        this.overlapRMSsimilarity = overlapRMSsimilarity;

        //  Calculating the threading energy according to the expanded sequence
        corpus.buildUngappedArray(seq.length);
        double[] threadingEnergies = calculateThreadingEnergy();

        // After the threading energy was calculated, the sequence can be cut to the essentials.
        this.seq = new int[fragEndInSeq - fragStartInSeq + 1];
        fragL = this.seq.length;
        for (int c = 0; c < fragL; c++)
            this.seq[c] = seq[fragStartInSeq + c];

        // And now building the lib
        similarToLib = new boolean[corpus.ungapped.length];
        for (int c = 0; c < corpus.ungapped.length; c++)
            similarToLib[c] = false;
        fragCoors = new double[3][3 * fragL];
        protCoors = new double[3][3 * fragL];
        System.out.println("Start building lib.");
        buildLib(threadingEnergies, fragStartInSeq, true);
        similarToLib = null;  // We don't need it anymore, and it takes quite a lot of memory

        // Where is the frag in the final sequence? i.e. not considering the overlap(s)
        if (manner == -1) {  // N term
            trueFragStart = overlap;
            trueFragEnd = fragL - 1;
        } else if (manner == 1) {  // C term
            trueFragStart = 0;
            trueFragEnd = fragL - 1 - overlap;
        } else if (manner == 0) {  // both sides
            trueFragStart = overlap;
            trueFragEnd = fragL - 1 - overlap;
        }
    }

    public LocalFragmentLib(Corpus corpus, int[] seq, int libSize, double rmsCutoff, int fragStartInSeq, int fragEndInSeq) {
        this.corpus = corpus;
        this.seq = seq;
        this.fragL = seq.length;
        this.libSize = libSize;
        this.rmsCutOff = rmsCutOff;

        //  Calculating the threading energy according to the expanded sequence
        corpus.buildUngappedArray(seq.length);
        double[] threadingEnergies = calculateThreadingEnergy();

        // After the threading energy was calculated, the sequence can be cut to the essentials.
        this.seq = new int[fragEndInSeq - fragStartInSeq + 1];
        fragL = this.seq.length;
        for (int c = 0; c < fragL; c++)
            this.seq[c] = seq[fragStartInSeq + c];

        similarToLib = new boolean[corpus.ungapped.length];
        for (int c = 0; c < corpus.ungapped.length; c++)
            similarToLib[c] = false;

        // And now building the lib
        similarToLib = new boolean[corpus.ungapped.length];
        for (int c = 0; c < corpus.ungapped.length; c++)
            similarToLib[c] = false;
        fragCoors = new double[3][3 * fragL];
        protCoors = new double[3][3 * fragL];
        System.out.println("Start building lib.");
        buildLib(threadingEnergies, fragStartInSeq, false);
        similarToLib = null;  // We don't need it anymore, and it takes quite a lot of memory
    }

    // This calculate the threading energy, according to the sequence and the the 'corpus.ungapped' array (built to the sequence length)

    private double[] calculateThreadingEnergy() {
        double[] energies = new double[corpus.ungapped.length];

        // Finding pre-Prolines in the sequence
        boolean[] pp = new boolean[seq.length];
        for (int e = 0; e < seq.length - 1; e++)
            if ((seq[e + 1] == 12) && (seq[e] != 5) && (seq[e] != 12))
                pp[e] = true;
            else
                pp[e] = false;
        pp[seq.length - 1] = false;

/*		// Random picking
		for (int c=0; c<corpus.ungapped.length ; c++) {
			energies[c] = Math.random();
		}
*/
/*		// If you want to calculate the energy according BLOSUM62, uncomment this loop and comment the next loop. 
        for (int c=0; c<corpus.ungapped.length ; c++) {
            energies[c] = 0.0;
            for (int e=0; e<seq.length ; e++) {
                energies[c] -= blosum62[seq[e]][corpus.resType[corpus.ungapped[c]+e]]; // Note the -= because we want to have energy not score
            }
        }
*/
        // Calculating the propensity threading energy
        for (int c = 0; c < corpus.ungapped.length; c++) {
            energies[c] = 0.0;
            for (int e = 0; e < seq.length; e++) {
                if (pp[e])
                    energies[c] += corpus.prePro[corpus.ungapped[c] + e];
                else
                    energies[c] += corpus.energies[corpus.ungapped[c] + e][0][seq[e]];
                energies[c] -= 0.2 * Corpus.blosum62[seq[e]][corpus.resType[corpus.ungapped[c] + e]]; // Note the -= because we want to have energy not score
            }
        }
        return energies;
    }

    private void buildLib(double[] energies, int fragStartInSeq, boolean filterStalks) {
        int[] tmpIndices = new int[libSize];
        double[] tmpEnergies = new double[libSize];

        for (int c = 0; c < libSize; c++) {
            tmpIndices[c] = -1;
            tmpEnergies[c] = 1e100;
        }

        // Start iterating for the lib
        for (int currentLibSize = 0; currentLibSize < libSize; currentLibSize++) {
            int lowInd = findLowestEnergy(tmpIndices, energies, filterStalks, fragStartInSeq);
            if (lowInd == -1) {
                System.out.println("WARNINIG!!! The corpus is too small to create such a large library. So Far created " + currentLibSize);
                libSize = currentLibSize;
                break;
            }
            tmpIndices[currentLibSize] = corpus.ungapped[lowInd] + fragStartInSeq;
            tmpEnergies[currentLibSize] = energies[lowInd];
            similarToLib[lowInd] = true;
            //System.out.println("Added to lib: " + currentLibSize);
        }

        // Updating the lib fields according to its final size
        libOrig = new int[libSize];
        libEnergy = new double[libSize];
        for (int c = 0; c < libSize; c++) {
            libOrig[c] = tmpIndices[c];
            libEnergy[c] = tmpEnergies[c];
        }
    }  // Of build lib


    public String printLib() {
        //if (calcRmsBetweenStruct(libOrig[0],corpus.ungapped[187]) < 1.0)
        //	return "Library discarded. Too close to a canonical helix.\n";
        String result = "";

        for (int i = 0; i < libSize; i++) {
            result += i + " " + libEnergy[i] + " " +
                    corpus.proteinNames[corpus.protInd[libOrig[i]]] + " " + corpus.resNum[libOrig[i]] + "\n";
        }

        return result;
    }

    public String printLib(int i) {
        //if (calcRmsBetweenStruct(libOrig[0],corpus.ungapped[187]) < 1.0)
        //	return "Library discarded. Too close to a canonical helix.\n";
        String result = "";

        result += i + " " + libEnergy[i] + " " +
                corpus.proteinNames[corpus.protInd[libOrig[i]]] + " " + corpus.resNum[libOrig[i]] + "\n";

        return result;
    }


    /**
     * This method checks if the fragment 'ind' (given in the ungapped indexing) is similar to a fragment already in the lib.
     * The array similarInLib is updated if a similarity is present.
     */
    private boolean checkSimilarity(int ind, int[] tmpIndices, int fragStartInSeq, boolean filterStalks) {
        if (similarToLib[ind])
            return true;
        for (int i = 0; tmpIndices[i] > -1; i++)
            if (!filterStalks) {
                if (corpus.calcRmsBetweenStruct(tmpIndices[i], corpus.ungapped[ind] + fragStartInSeq, fragL, 1, fragL) < rmsCutOff) {
                    similarToLib[ind] = true;
                    return true;
                }
            } else {
                if (corpus.calcRmsBetweenStruct(tmpIndices[i], corpus.ungapped[ind] + fragStartInSeq, fragL, manner, overlap) < rmsCutOff) {
                    similarToLib[ind] = true;
                    return true;
                }
            }
        return false;
    }

    /**
     * This method finds the lowest energy fragment that is not too similar to the rest of the lib *
     */
    private int findLowestEnergy(int[] tmpIndices, double[] energies, boolean filterStalks, int fragStartInSeq) {
        double lowestEnergy = 1e100;
        int lowestEnergyInd = -1;
        do {
            lowestEnergy = 1e100;
            lowestEnergyInd = -1;
            for (int i = 0; i < corpus.ungapped.length; i++) {
                if ((energies[i] < lowestEnergy) && !similarToLib[i]) {
                    if (!filterStalks || isFragCompatibleForConstruction(corpus.ungapped[i] + fragStartInSeq, overlapRMSsimilarity)) {
                        lowestEnergy = energies[i];
                        lowestEnergyInd = i;
                    } else
                        similarToLib[i] = true;
                }
            }
        } while ((lowestEnergyInd > -1) && checkSimilarity(lowestEnergyInd, tmpIndices, fragStartInSeq, filterStalks));
        return lowestEnergyInd;
    }


    /**
     * This method check whether the fragment 'indInLib' in the lib can overlap the given coordinates bellow a certain
     * RMS cut off. The fragment correspond to a protein fragment starting in residue 'residueFragStart'.
     * 'overlap' is the number of residues that are used for superposition. 'manner' means:
     * -1  - overlap of the N-terminus
     * 0 - two overlap regions in both termini (each of length overlap)
     * 1  - overlap of the C-terminus
     * <p/>
     * This method affects (and updates) the protCoors and fragCoors arrays.
     * 'prot' neednot be the protein from which the sequence for the library was taken.
     */
    public boolean isFragCompatible(int indInLib, Protein prot, int residueFragStart, int overlap, int manner, double overlapRMSsimilarity) {
        if (superimposeFrag(libOrig[indInLib], prot, residueFragStart, overlap, manner) < overlapRMSsimilarity)
            return true;
        else
            return false;
    }

    /**
     * The 'prot', 'overlap', 'manner' and 'residueFragStart' are from the class's fields.
     */
    public boolean isFragCompatible(int indInLib, double overlapRMSsimilarity) {
        if (prot == null)
            throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
        if (superimposeFrag(libOrig[indInLib], prot, residueFragStart, overlap, manner) < overlapRMSsimilarity)
            return true;
        else
            return false;
    }

    /**
     * This method is used only during the construction when 'libOrig' is not yet ready.
     */
    private boolean isFragCompatibleForConstruction(int indInCorpus, double overlapRMSsimilarity) {
        if (libOrig != null)
            throw new RuntimeException("\nERROR to the programer: you should call this method only if libOrig is still unavailable (i.e. ==null) because the constructor is running. \n");
        if (superimposeFrag(indInCorpus, prot, residueFragStart, overlap, manner) < overlapRMSsimilarity)
            return true;
        else
            return false;
    }


    /**
     * Superimposing fragment 'indInCorpus' in the corpus to a given portion in a protein.
     * The fragment correspond to a protein fragment starting in residue 'residueFragStart'.
     * 'overlap' is the number of residues that are used for superposition. 'manner' means:
     * -1  - overlap of the N-terminus
     * 0 - two overlap regions in both termini (each of length overlap)
     * 1  - overlap of the C-terminus
     * The RMS is calculated only isOn the overlap region.
     * <p/>
     * This method affects (and updates) the protCoors and fragCoors arrays.
     */
    public double superimposeFrag(int indInCorpus, Protein prot, int residueFragStart, int overlap, int manner) {
        // Calculating protCoor
        for (int c = 0; c < fragL; c++) {
            Atom atom = prot.residue(residueFragStart + c).amideN();
            protCoors[0][3 * c + 0] = atom.x();
            protCoors[1][3 * c + 0] = atom.y();
            protCoors[2][3 * c + 0] = atom.z();
            atom = prot.residue(residueFragStart + c).ca();
            protCoors[0][3 * c + 1] = atom.x();
            protCoors[1][3 * c + 1] = atom.y();
            protCoors[2][3 * c + 1] = atom.z();
            atom = prot.residue(residueFragStart + c).carbonylC();
            protCoors[0][3 * c + 2] = atom.x();
            protCoors[1][3 * c + 2] = atom.y();
            protCoors[2][3 * c + 2] = atom.z();
        }

        // Calculating fragCoor
        for (int c = 0; c < fragL; c++)
            for (int d = 0; d < 3; d++)
                for (int e = 0; e < 3; e++)
                    fragCoors[e][3 * c + d] = corpus.coors[indInCorpus + c][d][e];

        // Calculating the overlapping getAtoms according to manner
        int[] overlappingAtoms;
        if (manner == -1) {  // N term
            overlappingAtoms = new int[3 * overlap];
            for (int c = 0; c < 3 * overlap; c++)
                overlappingAtoms[c] = c;
        } else if (manner == 1) {  // C term
            overlappingAtoms = new int[3 * overlap];
            for (int c = 0; c < 3 * overlap; c++)
                overlappingAtoms[c] = 3 * fragL - 1 - c;
        } else if (manner == 0) {  // both sides
            overlappingAtoms = new int[2 * 3 * overlap];
            for (int c = 0; c < 3 * overlap; c++) {
                overlappingAtoms[c] = c;
                overlappingAtoms[3 * overlap + c] = 3 * fragL - 1 - c;
            }
        } else
            throw new RuntimeException("manner is not with valid value");
        try {
            double ret = Overlap.rmsPartialAltRMS(protCoors, fragCoors, overlappingAtoms);
            //System.out.println(ret);
            return ret;
        }
        catch (Exception e) {
            return 1e100;  // The overlap is not relaiable;
        }
    }


    /**
     * It is assumed that isFragCompatible was run before this method was called, so that fragCoors is updated.
     * The 'prot', 'overlap', 'manner' and 'residueFragStart' are from the class's fields.
     */
    public void insertFragToProt(DunbrackLib rotLib, int indInLib) {
        insertFragToProt(rotLib, indInLib, prot, residueFragStart, overlap, manner);
    }

    /**
     * It is assumed that isFragCompatible was run before this method was called, so that fragCoors is updated.
     * 'prot' neednot be the protein from which the sequence for the library was taken.
     */
    public void insertFragToProt(DunbrackLib rotLib, int indInLib, Protein prot, int residueFragStart, int overlap, int manner) {
        int startRes, endRes;
        if (manner == -1) {  // N term
            startRes = overlap;
            endRes = fragL - 1;
        } else if (manner == 1) {  // C term
            startRes = 0;
            endRes = fragL - 1 - overlap;
        } else if (manner == 0) {  // both sides
            startRes = overlap;
            endRes = fragL - 1 - overlap;
        } else
            throw new RuntimeException("manner is not with valid value");
        for (int c = startRes; c <= endRes; c++) {
            Atom atom = prot.residue(residueFragStart + c).amideN();
            atom.setXYZ(fragCoors[0][3 * c], fragCoors[1][3 * c], fragCoors[2][3 * c]);
            atom = prot.residue(residueFragStart + c).ca();
            atom.setXYZ(fragCoors[0][3 * c + 1], fragCoors[1][3 * c + 1], fragCoors[2][3 * c + 1]);
            atom = prot.residue(residueFragStart + c).carbonylC();
            atom.setXYZ(fragCoors[0][3 * c + 2], fragCoors[1][3 * c + 2], fragCoors[2][3 * c + 2]);
            ResidueBuilder.buildBackbone(prot.residue(residueFragStart + c),
                    corpus.torsions[libOrig[indInLib] + c][1],
                    corpus.torsions[libOrig[indInLib] + c][2]);
            // Into rot1
            double[] rot = null;
            if ((seq[c] != 0) && (seq[c] != 5))
                rot = rotLib.getRotamer(seq[c], corpus.torsions[libOrig[indInLib] + c][1],
                        corpus.torsions[libOrig[indInLib] + c][2], 0);
            ResidueBuilder.build(prot.residue(residueFragStart + c), seq[c], rot);
        }
    }


//	Getters...

    public int libSize() {
        return libSize;
    }

    public int libOrig(int ind) {
        return libOrig[ind];
    }

    public int trueFragStart() {
        if (prot == null)
            throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
        else
            return trueFragStart;
    }

    public int trueFragEnd() {
        if (prot == null)
            throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
        else
            return trueFragEnd;
    }

    public int residueFragStart() {
        if (prot == null)
            throw new RuntimeException("\nThis methods can be run only for instances that were initiated with constructor No. 1\n");
		else
			return residueFragStart; 
	}


} // Of LocalFragLib
