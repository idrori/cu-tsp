/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranCreator;
import meshi.geometry.rotamers.DunbrackLib;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;
import meshi.util.CommandList;
import meshi.util.UpdateableException;
import meshi.util.Utils;

import java.text.DecimalFormat;
import java.util.Vector;


/**
 * @author Nir 13/3/2008
 *         <p/>
 *         The purpose of this class is to correct models with impossible Ramachandran torsions, such as
 *         in the cases that occur after fragment assembly. The backbones of delinquent residues are replaced
 *         by other fragments in a manner that strives to keep the Ca's and Cb's as unpertubed as possible.
 *         <p/>
 *         The backbones of the immediately flanking residues of the rotten residue are also changed.
 *         <p/>
 *         The sidechains of the affected residues are placed in Rot1 (including the two flanking ones).
 */


public class CorrectRottenRamachandrans {

    // General parameters
    private final double CUTOFF_NON_GLY = 0.8; //1.3; Chen for debugging // A slightly less stringent value is 1.2 
    private final double CUTOFF_GLY = 0.9; //1.45;  Chen for debugging  // A slightly less stringent value is 1.35
    private final int MAX_LIB_SIZE_FOR_EACH_ROTTEN = 200;
    private final double DIVERSITY_RMS_IN_LIB = 0.1;
    private final double RMS_CUTOFF_IN_OVERLAP = 0.8;//0.4; Chen for debugging 
    private final double FACTOR_LARGE_CONTIG_RMS_CUTOFF_IN_OVERLAP = 0.2; // The final cutoff is: RMS_CUTOFF_IN_OVERLAP + (contig_length-1)*FACTOR_LARGE_CONTIG_RMS_CUTOFF_IN_OVERLAP


    // Local fields
    private Corpus corpus = null;
    private EnergyCreator[] creators = null;
    private DunbrackLib dunbrack = null;
    private CommandList commands = null;


    public CorrectRottenRamachandrans(CommandList commands, Corpus corpus) {
        this(commands, corpus, new DunbrackLib(commands, 1.0, 1));
    }

    /**
     * The Dunbrack lib is needed for Rot1 only.
     */
    public CorrectRottenRamachandrans(CommandList commands,
                                      Corpus corpus, DunbrackLib lib) {
        this.commands = commands;
        this.corpus = corpus;
        this.dunbrack = lib;  // For putting the sidechains into rot1 after altering
        creators = new EnergyCreator[1];
        creators[0] = new RamachandranCreator(1.0);
    }

    /**
     * This method returns an array with the residue numbers of the rotten-ramach residues.
     * NO changes to the coordinates are made.
     */

    public int[] detectRottenResidues(Protein prot) throws UpdateableException, EvaluationException {

        //FOR PRINTING when deducing CUTOFFs values:
        DecimalFormat fmt3 = new DecimalFormat("0.####");
        Vector<Integer> vec = new Vector<Integer>();

        TotalEnergy energy = new TotalEnergy(prot, creators, commands,"Generic total energy");
        try {
            energy.evaluate();
            energy.evaluateAtoms();
        } catch (UpdateableException ae) {
            throw new RuntimeException(ae);
        }
        double[] resRamachEnergies = new double[prot.residues().size()];
        for (int rr = 0; rr < resRamachEnergies.length; rr++) {
            if (!prot.residues().get(rr).dummy()) {
                resRamachEnergies[rr] = prot.residues().get(rr).ca().energy();
                if ((prot.residues().get(rr).type() == ResidueType.GLY) && (resRamachEnergies[rr] > CUTOFF_GLY)) {
                    vec.add(new Integer(prot.residues().get(rr).number()));
                    System.out.println(prot.residues().get(rr) + " " + fmt3.format(resRamachEnergies[rr]));
                }
                if ((prot.residues().get(rr).type() != ResidueType.GLY) && (resRamachEnergies[rr] > CUTOFF_NON_GLY)) {
                    vec.add(new Integer(prot.residues().get(rr).number()));
                    System.out.println(prot.residues().get(rr) + " " + fmt3.format(resRamachEnergies[rr]));
                }
            } else
                resRamachEnergies[rr] = Double.NEGATIVE_INFINITY;
        }

        // returning
        int[] result = new int[vec.size()];
        for (int c = 0; c < vec.size(); c++)
            result[c] = vec.get(c).intValue();

        return result;
    }

    /**
     * This correct the rotten-ramach residues as mentioned isOn top.
     * It returns an array with the residue numbers of the rotten-ramach residues.
     */
    public int[] detectAndCorrectRottenResidues(Protein prot) throws UpdateableException,EvaluationException{
        int[] rotten = detectRottenResidues(prot);

        Vector<Integer> contig = new Vector<Integer>();
        for (int pos = 0; pos < rotten.length; pos++) {
            contig.add(new Integer(rotten[pos]));
            correctFrag(prot, contig);
            contig = new Vector<Integer>();
        }

        return rotten;
    }

    private void correctFrag(Protein prot, Vector<Integer> contig) {
        LocalFragmentLib lib;
        double[][] savedCoordinates;
        int startRes = contig.get(0).intValue();  // Of the contig
        int endRes = contig.get(contig.size() - 1).intValue();  // Of the contig
        double bestRMS = Double.MAX_VALUE;
        int bestRMSfrag = -1;

        System.out.println("Correcting residues: ");
        for (int c = 0; c < contig.size(); c++)
            System.out.print(contig.get(c).intValue() + " ");

        // Can we correct also the flanking residues?
        // ------------------------------------------
        if ((prot.residue(startRes - 2) != null) && !prot.residue(startRes - 2).dummy() &&
                (prot.residue(endRes + 2) != null) && !prot.residue(endRes + 2).dummy()) { // Yes, we can
            savedCoordinates = saveCaCbCoordinates(prot, startRes - 1, endRes + 1);
            lib = new LocalFragmentLib(corpus, Utils.getSeqOfProt(prot, startRes - 2, endRes + 2), MAX_LIB_SIZE_FOR_EACH_ROTTEN * (contig.size() + 4),
                    DIVERSITY_RMS_IN_LIB, 0, endRes - startRes + 4,
                    prot, startRes - 2, 1, 0, RMS_CUTOFF_IN_OVERLAP + FACTOR_LARGE_CONTIG_RMS_CUTOFF_IN_OVERLAP * (contig.size() - 1));
            for (int fragNum = 0; fragNum < lib.libSize(); fragNum++) {
                lib.superimposeFrag(lib.libOrig(fragNum), prot, startRes - 2, 1, 0);
                lib.insertFragToProt(dunbrack, fragNum);
                if (calcRMS(prot, startRes - 1, endRes + 1, savedCoordinates) < bestRMS) {
                    bestRMSfrag = fragNum;
                    bestRMS = calcRMS(prot, startRes - 1, endRes + 1, savedCoordinates);
                }
            }
            lib.superimposeFrag(lib.libOrig(bestRMSfrag), prot, startRes - 2, 1, 0);
            lib.insertFragToProt(dunbrack, bestRMSfrag);
        } else {  // No, we cannot
            savedCoordinates = saveCaCbCoordinates(prot, startRes, endRes);
            lib = new LocalFragmentLib(corpus, Utils.getSeqOfProt(prot, startRes - 1, endRes + 1), MAX_LIB_SIZE_FOR_EACH_ROTTEN * (contig.size() + 2),
                    DIVERSITY_RMS_IN_LIB, 0, endRes - startRes + 2,
                    prot, startRes - 1, 1, 0, RMS_CUTOFF_IN_OVERLAP + FACTOR_LARGE_CONTIG_RMS_CUTOFF_IN_OVERLAP * (contig.size() - 1));
            for (int fragNum = 0; fragNum < lib.libSize(); fragNum++) {
                lib.superimposeFrag(lib.libOrig(fragNum), prot, startRes - 1, 1, 0);
                lib.insertFragToProt(dunbrack, fragNum);
                if (calcRMS(prot, startRes, endRes, savedCoordinates) < bestRMS) {
                    bestRMSfrag = fragNum;
                    bestRMS = calcRMS(prot, startRes, endRes, savedCoordinates);
                }
            }
            lib.superimposeFrag(lib.libOrig(bestRMSfrag), prot, startRes - 1, 1, 0);
            lib.insertFragToProt(dunbrack, bestRMSfrag);
        }

        System.out.println("     Best RMS: " + bestRMS);
    }


    private double calcRMS(Protein prot, int start, int end, double[][] refCoors) {
        double totRms = 0.0;
        Atom atom;
        int c = 0;
        for (int res = start; res <= end; res++) {
            atom = prot.residue(res).ca();
            totRms += (atom.x() - refCoors[0][c]) *
                    (atom.x() - refCoors[0][c]) +
                    (atom.y() - refCoors[1][c]) *
                            (atom.y() - refCoors[1][c]) +
                    (atom.z() - refCoors[2][c]) *
                            (atom.z() - refCoors[2][c]);
            c++;
            if (prot.residue(res).cb() != null) {
                atom = prot.residue(res).cb();
                totRms += (atom.x() - refCoors[0][c]) *
                        (atom.x() - refCoors[0][c]) +
                        (atom.y() - refCoors[1][c]) *
                                (atom.y() - refCoors[1][c]) +
                        (atom.z() - refCoors[2][c]) *
                                (atom.z() - refCoors[2][c]);
                c++;
            }
        }
        return Math.sqrt(totRms / c);
    }


    private double[][] saveCaCbCoordinates(Protein prot, int resStart, int resEnd) {
        double[][] savedLoopCoordinates = new double[3][(resEnd - resStart + 1) * 2];
        Atom atom;
        int counter = 0;
        for (int c = resStart; c <= resEnd; c++) {
            atom = prot.residue(c).ca();
            savedLoopCoordinates[0][counter] = atom.x();
            savedLoopCoordinates[1][counter] = atom.y();
            savedLoopCoordinates[2][counter] = atom.z();
            counter++;
            if (prot.residue(c).cb() != null) {
                atom = prot.residue(c).cb();
                savedLoopCoordinates[0][counter] = atom.x();
                savedLoopCoordinates[1][counter] = atom.y();
                savedLoopCoordinates[2][counter] = atom.z();
                counter++;
            }
        }
        return savedLoopCoordinates;
    }


}  // of CorrectRottenRamachandrans
