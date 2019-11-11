/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.rg.filters.HeavyAtoms;
import meshi.util.overlap.*;
import meshi.util.filters.*;
import meshi.sequences.*;
import meshi.geometry.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.ca.*;

import java.util.*;
import java.lang.*;


/**
 * A (hopefully) comfortable handle to the complexity of the Kabsch structural overlap algorithm[1,2].
 * The algorithm itself is implemented in meshi.overlap.Overlap.
 * 1. "A solution for the best rotation to relate two sets of vectors". By Wolfgang Kabsch, (1976) Acta Cryst. A 32: 922-923<BR>
 * 2. "A discussion of the solution for the best rotation to relate two sets<BR>
 * of vectors". By Wolfgang Kabsch, Acta Cryst. (1978). A34, 827-828.
 * This class was rather drastically modified by Chen in version 1.45.
 */
public class Rms implements KeyWords {

    private static boolean debug = true;
    private boolean alive = true;
    //data members
    private double[][] rotateMatrix;
    private double rms;
    private Coordinates centerOfMass0 = new Coordinates();
    private Coordinates centerOfMass1 = new Coordinates();
    private static double[] zeros = {0.0, 0.0, 0.0, 0.0, 0.0};
    public static enum RmsType {CA, ALL_ATOMS, HEAVY_ATOMS}
    public static final Filter caFilter = new CaFilter();
    public static final Filter kolDicfin = new KolDichfin();
    public static final Filter heavyFilter = new HeavyAtoms();
    //constructors

    public Rms(ResidueAlignment residueAlignment) {
        this(new AtomAlignment(ResidueAlignment.toAlignmentArray(residueAlignment), new CaFilter()));
    }



    public Rms(AtomAlignment atomAlignment) {
        if (atomAlignment.hasGaps()) throw new RuntimeException("Cannot calculate RMS for AtomAlignment with gaps");
        int size = atomAlignment.size();
        double[][] coor0 = new double[3][size];
        double[][] coor1 = new double[3][size];
        String comment0 = atomAlignment.get(0).cell(0).comment;
        String comment1 = atomAlignment.get(0).cell(1).comment;
        double centerOfMassX1 = 0, centerOfMassY1 = 0, centerOfMassZ1 = 0;
        double centerOfMassX0 = 0, centerOfMassY0 = 0, centerOfMassZ0 = 0;
        try {
            for (int i = 0; i < size; i++) {
                Atom atom0 = atomAlignment.atomAt(i, 0);
                Atom atom1 = atomAlignment.atomAt(i, 1);
                coor0[0][i] = atom0.x();
                coor0[1][i] = atom0.y();
                coor0[2][i] = atom0.z();
                coor1[0][i] = atom1.x();
                coor1[1][i] = atom1.y();
                coor1[2][i] = atom1.z();

                centerOfMassX0 += atom0.x();
                centerOfMassY0 += atom0.y();
                centerOfMassZ0 += atom0.z();
                centerOfMassX1 += atom1.x();
                centerOfMassY1 += atom1.y();
                centerOfMassZ1 += atom1.z();
            }
            centerOfMassX0 /= size;
            centerOfMassY0 /= size;
            centerOfMassZ0 /= size;
            centerOfMassX1 /= size;
            centerOfMassY1 /= size;
            centerOfMassZ1 /= size;

            centerOfMass0.setXYZ(centerOfMassX0, centerOfMassY0, centerOfMassZ0);
            centerOfMass1.setXYZ(centerOfMassX1, centerOfMassY1, centerOfMassZ1);

            Overlap overlap = new Overlap(coor0, coor1, size, comment0, comment1);
            rms = overlap.rms();
            rotateMatrix = overlap.rotationMatrix();
            if (Double.isNaN(rotateMatrix[0][0]))
                throw new RuntimeException("This is weird.");
        }
        catch (Exception e) {
            throw new MeshiException("Rms Error:\n" +
                    "comparing " + comment0 + "\n" +
                    "with\n" + comment1 + "\n" +
                    "reproted problem" +
                    e);
        }
    }


    private Coordinates centerOfMass0() {
        return centerOfMass0;
    }

    private Coordinates centerOfMass1() {
        return centerOfMass1;
    }

    public static double superimpose(Protein protein0, Protein protein1, ResidueAlignment residueAlignment) {

        ResidueAlignment[] residueAlignments = ResidueAlignment.toAlignmentArray(residueAlignment);
        AtomAlignment atomAlignment = new AtomAlignment(residueAlignments, new CaFilter());
        Rms rms = new Rms(atomAlignment);
        Coordinates centerOfMass0 = rms.centerOfMass0();
        Coordinates centerOfMass1 = rms.centerOfMass1();
        for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            atom.addToX(0 - centerOfMass1.x());
            atom.addToY(0 - centerOfMass1.y());
            atom.addToZ(0 - centerOfMass1.z());
        }
        double[][] rotateMatrix = rms.getMatrix();
        for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            double x = atom.x();
            double y = atom.y();
            double z = atom.z();
            atom.setXYZ(rotateMatrix[0][0] * x + rotateMatrix[0][1] * y + rotateMatrix[0][2] * z,
                    rotateMatrix[1][0] * x + rotateMatrix[1][1] * y + rotateMatrix[1][2] * z,
                    rotateMatrix[2][0] * x + rotateMatrix[2][1] * y + rotateMatrix[2][2] * z);
        }
        for (Iterator atoms = protein1.atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            atom.addToX(centerOfMass0.x());
            atom.addToY(centerOfMass0.y());
            atom.addToZ(centerOfMass0.z());
        }
        return rms.getRms();
    }


    //methods

    public double[][] getMatrix() {
        if (alive) return rotateMatrix;
        else return null;
    }


    public double getRms() {
        if (alive) return rms;
        else return -1;
    }

    public String toString() {
        return (new Double(getRms())).toString();
    }

    public Coordinates getCenterOfMass0() {
        return centerOfMass0;
    }

    public Coordinates getCenterOfMass1() {
        return centerOfMass1;
    }

    //------------------------------------------------------------------ utility methods --------------------------------------------------------------

    public static final double[] NO_GDT = {-1, -1, -1, -1, -1};

    public static double rms(Protein protein0, CommandList commands,
                             ResidueAlignmentMethod residueAlignmentMethod) throws AlignmentException{
        return rms(protein0, commands, new KolDichfin(), residueAlignmentMethod); // KolDichfin - a filter that accept all
    }

    public static double rms(Protein protein0,
                             CommandList commands,
                             Filter filter,
                             ResidueAlignmentMethod residueAlignmentMethod)throws AlignmentException{
        AtomAlignment atomAlignment = getAtomAlignment(protein0, commands, filter, residueAlignmentMethod);
        if (atomAlignment == null) return -1;
        Rms rms = new Rms(atomAlignment);
        return rms.getRms();
    }

    public static double rms(Protein protein0, Protein protein1,
                             ResidueAlignmentMethod residueAlignmentMethod)throws AlignmentException{
        return rms(protein0,protein1,residueAlignmentMethod,RmsType.CA);
    }
    public static double rmsAll(Protein protein0, Protein protein1,
                                ResidueAlignmentMethod residueAlignmentMethod) throws AlignmentException{
            return rms(protein0,protein1,residueAlignmentMethod,RmsType.ALL_ATOMS);
        }
    public static double rmsHeavy(Protein protein0, Protein protein1,
                                  ResidueAlignmentMethod residueAlignmentMethod) throws AlignmentException{
                return rms(protein0,protein1,residueAlignmentMethod,RmsType.HEAVY_ATOMS);
            }

    public static double rms(Protein protein0, Protein protein1,
                             ResidueAlignmentMethod residueAlignmentMethod,RmsType type) throws AlignmentException{
        ResidueAlignment alignment = new ResidueAlignment(protein0.chains(), protein0.name(), protein1.chains(), protein1.name(), residueAlignmentMethod);
        if (alignment.size() < 5) {
            String message = "Warrning: Cannot calculate rms for "
                    + protein0.name() + " and " + protein1.name()
                    + " using "+residueAlignmentMethod + "\n"
                    + "alignment.size() < 5 " + alignment.size()+" returning RMS of -1";
            System.out.println(message);
            System.err.println(message);
        }
        return rms(alignment,type);
    }

    public static double rms(Protein protein0,
                             Protein protein1,
                             ResidueAlignmentMethod residueAlignmentMethod,
                             Filter filter) throws AlignmentException{
        ResidueAlignment alignment = new ResidueAlignment(protein0.chains(), protein0.name(), protein1.chains(), protein1.name(), residueAlignmentMethod);
        ResidueAlignment filteredAlignment = new ResidueAlignment();
        for (ResidueAlignmentColumn column : alignment) {
            if (filter.accept(column)) filteredAlignment.add(column);
        }
        return rms(filteredAlignment,RmsType.CA);
    }





    //modified constructors Added by Tommer 1.9.14:

    public static double[] gdt(Protein protein0, Protein protein1){
        return  gdt(protein0, protein1, GDTcalculator.Type.TS);
    }



    public static double[] gdt(Protein protein0, Protein protein1, //Added By Tommer 1.9.14
			                   GDTcalculator.Type type){
    	//long startTime = System.currentTimeMillis();
        ResidueAlignment residueAlignment = null;
        try {
            ChainList chains0 = protein0.chains();
            ChainList chains1 = protein1.chains();
            String    name0   = protein0.name();
            String    name1   = protein1.name();
            residueAlignment = new ResidueAlignment(chains0, name0, chains1, name1, ResidueAlignmentMethod.IDENTITY);//modified
        }
        catch (AlignmentException ex) {
            Utils.throwException("Static function Rms.gdt",ex,"Failed to align "+protein0+" and "+protein1);
        }
		//long estimatedTime = System.currentTimeMillis() - startTime;
		//System.out.println("\nestimated time: "+estimatedTime);

				if (residueAlignment.size() < 5) {
				    String message = "Warrning: Cannot calculate GDT_TS for "
                            + protein0.name() + " and " + protein1.name()
                            + " using identity alignment + \n"
                            + "residueAlignment.size() < 5 " + residueAlignment.size()+" returning GDT of 0";
                    System.out.println(message);
                    System.err.println(message);
                    return new double[5];
                }
        //residueAlignment.print();
		return gdt(residueAlignment, protein0.atoms().CAFilter().size(), type);
	}

  //END OF modified constructors Added by Tommer 1.9.14

    public static double[] gdt(Protein protein0, Protein protein1, Filter filter)throws AlignmentException {
        return gdt(protein0, protein1,filter,GDTcalculator.Type.TS);
    }
    public static double[] gdt(Protein protein0,
                               Protein protein1,
                               Filter filter, GDTcalculator.Type type) throws AlignmentException{
            ResidueAlignment residueAlignment = new ResidueAlignment(protein0.chains(), protein0.name(),
                protein1.chains(), protein1.name(),
                ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
        ResidueAlignment filteredAlignment = new ResidueAlignment();
        for (ResidueAlignmentColumn column : residueAlignment) {
            if (filter.accept(column)) filteredAlignment.add(column);
        }
        return gdt(filteredAlignment, protein0.atoms().CAFilter().size(), type);
    }


    public static double[] gdt(ResidueAlignment residueAlignment, int refLength) {
        return gdt(residueAlignment,refLength, GDTcalculator.Type.TS);
    }
        public static double[] gdt(ResidueAlignment residueAlignment, int refLength, GDTcalculator.Type type) {
            return gdt(residueAlignment, new KolDichfin(), refLength, type);
    }

    public static double rms(ResidueAlignment residueAlignment, RmsType type) {
        Filter filter;
        switch (type){
            case CA: filter = caFilter; break;
            case HEAVY_ATOMS: filter = heavyFilter; break;
            case ALL_ATOMS: filter = kolDicfin; break;
            default: throw new RuntimeException("Weird RmsType "+type);
        }
        ResidueAlignment[] residueAlignments = ResidueAlignment.toAlignmentArray(residueAlignment);
        AtomAlignment atomAlignment = new AtomAlignment(residueAlignments, filter);
        if (atomAlignment == null) return -1;
        Rms rms = new Rms(atomAlignment);
        return rms.getRms();
    }

    public static double[] gdt(Protein protein0, CommandList commands) throws AlignmentException{
        return gdt(protein0, commands, new KolDichfin()); // KolDichfin - a filter that accept all
    }

    public static double[] gdt(Protein model, CommandList commands, Filter filter) throws AlignmentException {
        //	try {
        CommandList alignmentCommands = commands.firstWordFilter(SUPERIMPOSE);
        String referenceFileName = alignmentCommands.secondWord(REFERENCE).thirdWord();
        if (referenceFileName.equals(NONE.key)) {
            System.out.println("No reference structure");
            return NO_GDT;
        }
        Protein reference = Utils.getProtein(commands, SUPERIMPOSE, REFERENCE, CaResidue.creator, Utils.defaultExceptionHandler);
        int refLength = reference.atoms().CAFilter().size();
        ResidueAlignment residueAlignment = new ResidueAlignment(reference.chains(), reference.name(), model.chains(), model.name(), ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
        return gdt(residueAlignment, filter, refLength);
    }

    public static double[] gdt(ResidueAlignment residueAlignment, Filter filter, int refLength) {
         return gdt(residueAlignment,filter,refLength, GDTcalculator.Type.TS);
    }
        public static double[] gdt(ResidueAlignment residueAlignment, Filter filter, int refLength, GDTcalculator.Type type) {
           AtomList modelAtoms;
            ResidueAlignment[] residueAlignments = ResidueAlignment.toAlignmentArray(residueAlignment);
            AtomAlignment atomAlignment = new AtomAlignment(residueAlignments, filter);
        if (atomAlignment == null) return NO_GDT;
        AtomList referenceAtoms = atomAlignment.atomList(0);
        try {
            modelAtoms = atomAlignment.atomList(1);
        }
        catch (Exception ex) {
            System.out.println("Problem while getting an AtomList from the atomAlignment");
            atomAlignment.print();
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
        AtomList newReferenceAtoms = Utils.duplicateInAnewMolecularSystem(referenceAtoms,new MolecularSystem());
        AtomList newModelAtoms = Utils.duplicateInAnewMolecularSystem(modelAtoms, new MolecularSystem());

        return GDTcalculator.gdt(newReferenceAtoms, newModelAtoms, refLength,type);
    }


    private static AtomAlignment getAtomAlignment(Protein protein0,
                                                  CommandList commands,
                                                  Filter filter,
                                                  ResidueAlignmentMethod
                                                          residueAlignmentMethod) throws AlignmentException {
        CommandList alignmentCommands = commands.firstWordFilter(SUPERIMPOSE);
        String referenceFileName = alignmentCommands.secondWord(REFERENCE).thirdWord();
        if (referenceFileName.equals(NONE.key)) {
            System.out.println("No reference structure");
            return null;
        } else {
            Protein protein1 = Utils.getProtein(commands, SUPERIMPOSE, REFERENCE, CaResidue.creator, Utils.defaultExceptionHandler);
            if (protein0.chains().size() != protein1.chains().size())
                throw new RuntimeException("Cannot compare " + protein0 + " with " + protein0.chains().size() + " chains and " + protein1 + " with " + protein1.chains().size() + " chanins.");
            if (debug) System.out.println("protein1 = " + protein1);
            String mode = alignmentCommands.secondWord(RMS_MODE).thirdWord();
            if (!mode.equals(ALL_CA.key))
                throw new RuntimeException("Currently the only implemented mode of " + SUPERIMPOSE + "\n" +
                        "is " + ALL_CA);
            MeshiSequence[] sequences0 = protein0.sequences();
            MeshiSequence[] sequences1 = protein1.sequences();
            ResidueAlignment[] residueAlignments = new ResidueAlignment[protein0.chains().size()];
            for (int iChain = 0; iChain < protein0.chains().size(); iChain++) {
                Utils.print("Generating ResidueAlignment of \n" + sequences0[iChain] + "\n" + sequences1[iChain]);
                residueAlignments[iChain] = new ResidueAlignment(protein1.chains().get(iChain), protein1.name(), protein0.chains().get(iChain), protein0.name(), residueAlignmentMethod);
                Utils.print(residueAlignments[iChain]);
            }
            return new AtomAlignment(residueAlignments, filter);
        }
    }
}
