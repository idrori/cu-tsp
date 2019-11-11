/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.applications.prediction.OriginalAtom;
import meshi.energy.compatebility.ResidueSsPrediction;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSumma;
import meshi.energy.simpleEnergyTerms.GoodBonds;
import meshi.energy.simpleEnergyTerms.ParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergy;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandran.RamachandranEnergyElement;
import meshi.sequences.aligner.IdentityMatrix;
import meshi.util.file.*;
import meshi.util.filters.*;
import meshi.util.dssp.*;
import meshi.util.string.*;
import meshi.geometry.*;
import meshi.PDB.*;
import meshi.energy.*;
import meshi.optimizers.*;
import meshi.parameters.*;
import meshi.sequences.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.ca.*;
import meshi.molecularElements.extendedAtoms.*;
import meshi.energy.simpleEnergyTerms.angle.*;
import meshi.energy.simpleEnergyTerms.plane.*;
import meshi.energy.simpleEnergyTerms.outOfPlane.*;
import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.*;
//import meshi.applications.prediction.*;
import java.util.*;
import java.io.*;
import javax.media.j3d.Transform3D;
import javax.vecmath.*;

/**
 * Where we store useful static methods that do not make sense anywhere else.
 */

public class Utils implements MeshiPotential, KeyWords {
    public static final ExceptionHandler defaultExceptionHandler = new DefaultExceptionHandler();
    public static final ExceptionHandler doNothing = new DoNothingExceptionHandler();

    //----------------------------------- Commands -------------------------------------------------

    /**
     * A standard method to create the commands list.
     */

    public static CommandList init(String[] args, String name) {
        return init(args, 1, ("Usage: java -XmxNNNm " + name + "  <commands file name>\n\n" +
                "NNN is the size of the expected memory requirement in MegaBytes."));
    }

    public static CommandList init(String[] args, int numberOfParameters, String errorMessage) {
        if ((numberOfParameters > 0) && (args.length != numberOfParameters)) throw new RuntimeException(errorMessage);
        CommandList commands = new CommandList(args,
                new CommandsException(errorMessage));
        if (commands.keyExists(SEED))
            MeshiProgram.initRandom(commands.firstWord(SEED).secondWordInt());
        return commands;
    }

    public static CommandList init(String[] args, int numberOfParameters, int seed, String errorMessage) {
        if ((numberOfParameters > 0) && (args.length != numberOfParameters)) throw new RuntimeException(errorMessage);
        CommandList commands = new CommandList(args,
                new CommandsException(errorMessage));
        MeshiProgram.initRandom(seed);
        return commands;
    }


    //----------------------------------- Alignments ----------------------------------------------

    /**
     * Gets a residue alignment as a parameters and returns a new alignment that includes only residues with original getAtoms.
     */
    public static ResidueAlignment getOriginalAtomsAlignment(ResidueAlignment alignment, OriginalAtoms originalAtoms,
                                                             int row) {
        ResidueAlignment out = new ResidueAlignment();
        for (Iterator columns = alignment.iterator(); columns.hasNext();) {
            ResidueAlignmentColumn column = (ResidueAlignmentColumn) columns.next();
            Residue residue = (Residue) column.cell(row).obj;
            if (!residue.dummy()) {
                Atom atom = residue.ca();
                if (originalAtoms.accept(atom)) out.add(column);
            }
        }
        return out;
    }


    //----------------------------------- Updateable -------------------------------------------
     public static void checkUpdate(int resourceNumberOfUpdates,
                                    int requestedNumberOfUpdates,
                                    Updateable resource) {
        if (requestedNumberOfUpdates != resourceNumberOfUpdates + 1)
            throw new RuntimeException("Something weird with "+resource+".update(int numberOfUpdates)\n" +
                    "requested numberOfUpdates = " + requestedNumberOfUpdates + "  "+
                    resource+" numberOfUpdates = " + resourceNumberOfUpdates);
    }
    //----------------------------------- DistanceMatrix -------------------------------------------

    /**
     * A standard method to build distance matrices.

    public static DistanceMatrix getDistanceMatrix(AtomList atomList, CommandList commands) {
        double rmax = commands.firstWord(R_MAX).secondWordDouble();
        double buffer = commands.firstWord(BUFFER_SIZE).secondWordDouble();
        double edge = commands.firstWord(GRID_EDGE).secondWordDouble();
        return new DistanceMatrix(atomList.molecularSystem(), rmax, buffer, edge, DistanceMatrix.DEFAULT_BONDED_LIST_DEPTH);
    }
     */
    //----------------------------------- energy -------------------------------------------------

    public static double[] getEnergyValues(TotalEnergy totalEnergy) {
        double[] out = new double[totalEnergy.energyTerms().size()];
        boolean[] on = new boolean[totalEnergy.energyTerms().size()];
        int i = 0;
        for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
            AbstractEnergy term = (AbstractEnergy) terms.next();
            on[i] = term.isOn();
            i++;
        }
        i = 0;
        for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
            AbstractEnergy term = (AbstractEnergy) terms.next();
            totalEnergy.off();
            term.on();
            try {
                out[i] = totalEnergy.evaluate();
            } catch (Exception ex) {
                System.out.println(" energy.evaluate failed die to " + ex);
                ex.printStackTrace();
                throw new RuntimeException("Quiting");
            }
            i++;
        }
        i = 0;
        totalEnergy.off();
        for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
            AbstractEnergy term = (AbstractEnergy) terms.next();
            if (on[i]) term.on();
            i++;
        }
        return out;
    }

    public static String[] getEnergyTermsNames(TotalEnergy totalEnergy) {
        String[] out = new String[totalEnergy.energyTerms().size()];
        int i = 0;
        for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
            AbstractEnergy term = (AbstractEnergy) terms.next();
            out[i] = term.comment();
            i++;
        }
        return out;
    }

    public static TotalEnergy relax(Protein protein,
                                    EnergyCreator[] energyCreators,
                                    CommandList commands)throws UpdateableException, EvaluationException,AlignmentException{
        return relax(protein.atoms(),protein,energyCreators,commands);
    }

        public static TotalEnergy relax(AtomList atoms,
                                    Protein protein,
                                    EnergyCreator[] energyCreators,
                                    CommandList commands)throws UpdateableException, EvaluationException,AlignmentException{
        return relax(atoms, protein, energyCreators, DistanceMatrix.DistanceMatrixType.STANDARD,
                     commands,RELAX);
    }

    public static TotalEnergy relax(AtomList atoms,
                                    Protein protein,
                                    EnergyCreator[] energyCreators,DistanceMatrix.DistanceMatrixType distanceMatrixType,
                                    CommandList commands)throws UpdateableException, EvaluationException,AlignmentException{
        return relax(atoms, protein, energyCreators, distanceMatrixType,commands, RELAX);
    }

    public static TotalEnergy oldRelax(AtomList atoms, Protein protein, EnergyCreator[] energyCreators, CommandList commands)throws UpdateableException,EvaluationException{
        return oldRelax(atoms, protein, energyCreators, commands, RELAX);
    }

    public static TotalEnergy oldRelax(AtomList atoms,
                                       Protein protein, EnergyCreator[] energyCreators,
                                       CommandList commands, Key key) throws UpdateableException,EvaluationException{
        // chen 3/3/12

        TotalEnergy energy = new TotalEnergy(protein, energyCreators, commands,"generic energy (oldRelax) ");
        energy.evaluate();
        Optimizer optimizer = Utils.getLBFGS(energy, commands, key);
        try {
            optimizer.run();
        }
        catch (Exception ex) {
            try {
                energy.test();
            } catch (Exception ex1) {
                System.out.println(" energy.test failed die to " + ex);
                ex.printStackTrace();
                throw new RuntimeException("Quiting");
            }
            throw new RuntimeException(ex);
        }
        return energy;
    }

    public static TotalEnergy relax(AtomList atoms, Protein protein,
                                    EnergyCreator[] energyCreators, DistanceMatrix.DistanceMatrixType distanceMatrixType,
                                    CommandList commands, Key key) throws UpdateableException, EvaluationException,AlignmentException{

        TotalEnergy energy = new TotalEnergy(protein, energyCreators, distanceMatrixType, commands,"generic energy (relax)");
        Optimizer optimizer = Utils.getLBFGS(energy, commands, key);
        AtomicPairwisePMFSumma summaTerm = (AtomicPairwisePMFSumma) energy.getEnergyTerm(new AtomicPairwisePMFSumma());
        if ((summaTerm != null) && (summaTerm.weight() > 1)) summaTerm.setWeight(1);
        try {
            optimizer.run();
        }
        catch (OptimizerException ex) {
            if (verbose) {
                try {
                    energy.test();
                } catch (Exception ex1) {
                    System.out.println(" energy.evaluate failed die to " + ex);
                    ex.printStackTrace();
                    throw new RuntimeException("Quiting");
                }
                throw new RuntimeException(ex);
            }
            else {System.out.println(" energy.evaluate failed die to " + ex);
                    ex.printStackTrace();
                    throw new RuntimeException("Quiting");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        energy.evaluateAtoms();
        return energy;
    }

    //----------------------------------- io and files ----------------------------------------------------

    public static StringList getStructureNames(CommandList commands) {
        String fileName = commands.firstWord(STRUCTURE_NAMES).secondWord();
        return new StringList(new MeshiLineReader(fileName));
    }

    public static MeshiWriter getOutputFile(CommandList commands) {
        String fileName = commands.firstWord(OUTPUT_FILE_NAME).secondWord();
        try {
            return new MeshiWriter(fileName);
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
    }

    /**
     * Opens a file for writing (warped by MeshiWriter object) with a name and location specified
     * by the commands and a .<seed>.<runNumber>.pdb suffix.
     */
    public static MeshiWriter newPdbWriter(CommandList commands, Key pathKey, Key nameKey, int runNumber, int templateNumber) {
        String templateNumberString;
            try {
                if (templateNumber < 0) {
                    templateNumberString = "";
                }
                else templateNumberString = "."+templateNumber;
                String suffix = "." + MeshiProgram.seed() + "." + runNumber +templateNumberString+ ".pdb";                                                                                                      //MOTY
                return new MeshiWriter(commands, pathKey, nameKey, suffix);
            }
            catch (Exception ex) {
                throw new RuntimeException("A problem in opening output file with key words " +
                        pathKey + " and " + nameKey + "\n" + ex);
            }
        }

    public static MeshiWriter newPdbWriter(CommandList commands, Key pathKey, Key nameKey, int runNumber) {
        return newPdbWriter(commands,pathKey,nameKey,runNumber,-1);
    }

    public  static ArrayList<File>   getDirectories(String prefix) {
        File thisDir = new File(".");
        File[] files = thisDir.listFiles();
        ArrayList<File> directories = new ArrayList<File>();
        for (File file : files) {
            if (file.getName().startsWith(prefix) && file.isDirectory())
                  directories.add(file);
        }
        return directories;
    }

    //----------------------------------- analysis ------------------------------------------------

    public static void RmsGdtEnergy(Object output,
                                    Chain model,
                                    Chain reference,
                                    TotalEnergy totalEnergy,
                                    String tag) throws EvaluationException,AlignmentException{
        ResidueAlignment modelRefeenceAlignment = null;
        int referenceSize = -1;
        if (reference != null) {
            modelRefeenceAlignment = new ResidueAlignment(model,"model", reference,"reference",ResidueAlignmentMethod.BY_RESIDUE_NUMBER);
            referenceSize = reference.size();
        }
        if (output instanceof StringList)
            RmsGdtEnergy((StringList) output, modelRefeenceAlignment, totalEnergy, tag, referenceSize, true);
        else if (output instanceof MeshiWriter)
            RmsGdtEnergy((MeshiWriter) output, modelRefeenceAlignment, totalEnergy, tag, referenceSize, true);
        else throw new RuntimeException("Wrong first parameter for RmsGdtEnergy " + output);
    }

    public static void RmsGdtEnergy(Object output, ResidueSequence model,
                                    ResidueSequence reference, TotalEnergy totalEnergy,
                                    String tag) throws EvaluationException{
        ResidueAlignment modelRefeenceAlignment = null;
        int referenceSize = -1;
        if (reference != null) {
            modelRefeenceAlignment = new ResidueAlignment(model, reference, new IdentityMatrix());
            referenceSize = reference.size();
        }
        if (output instanceof StringList)
            RmsGdtEnergy((StringList) output, modelRefeenceAlignment, totalEnergy, tag, referenceSize, true);
        else if (output instanceof MeshiWriter)
            RmsGdtEnergy((MeshiWriter) output, modelRefeenceAlignment, totalEnergy, tag, referenceSize, true);
        else throw new RuntimeException("Wrong first parameter for RmsGdtEnergy " + output);
    }




    public static void RmsGdtEnergy(MeshiWriter output,
                                    ResidueAlignment modelReferenceAlignment,
                                    TotalEnergy totalEnergy, String tag,
                                    int refLength, boolean headerFlag) throws EvaluationException{
        StringList stringList = new StringList();
        RmsGdtEnergy(stringList, modelReferenceAlignment, totalEnergy, tag, refLength, headerFlag);
        stringList.print(output);
        output.flush();
    }

    public static void RmsGdtEnergy(StringList output,
                                    ResidueAlignment modelReferenceAlignment,
                                    TotalEnergy totalEnergy, String tag,
                                    int refLength, boolean headerFlag) throws EvaluationException{
        double rms = -1;
        double rmsOrig = -1;
        double[] gdt = {-1.0, -1.0, -1.0, -1.0, -1.0};
        double[] gdtOrig = {-1.0, -1.0, -1.0, -1.0, -1.0};
        String header = "H_" + tag;
        String values = "V_" + tag;
        if (modelReferenceAlignment != null) {
            ResidueAlignment origRefeenceAlignment =  modelReferenceAlignment.filter(OriginalAtom.filter);
            rms = Rms.rms(modelReferenceAlignment,Rms.RmsType.CA);
            if (origRefeenceAlignment.size() > 3)
                rmsOrig = Rms.rms(origRefeenceAlignment,Rms.RmsType.CA);
            gdt = Rms.gdt(modelReferenceAlignment, refLength);
            if (origRefeenceAlignment.size() > 3)
                gdtOrig = Rms.gdt(origRefeenceAlignment, refLength);
            if (gdtOrig[0] > 0) {
                header += String.format("%6s %6s %5s %5s  %5s  %5s  %5s  %5s ", "RMS_O", "RMS", "GDT_O", "GDT_TS", "GDT1", "GDT`2", "GDT3", "GDT4");
                values += String.format(" %6.3f %6.3f %5.3f %5.3f  %5.3f  %5.3f  %5.3f  %5.3f ",
                        rmsOrig, rms, gdtOrig[0], gdt[0], gdt[1], gdt[2], gdt[3], gdt[4]);
            } else {
                header += String.format("%6s %5s %5s  %5s  %5s  %5s ", "RMS", "GDT_TS", "GDT1", "GDT2", "GDT3", "GDT4");
                values += String.format("%6.3f %5.3f  %5.3f  %5.3f  %5.3f  %5.3f ",
                        rms, gdt[0], gdt[1], gdt[2], gdt[3], gdt[4]);
            }
        }
        double energy;
        if (totalEnergy != null) {
            try {
                energy = totalEnergy.evaluate();
            } catch (Exception ae) {
                throw new RuntimeException("RmsGdtEnergy failed due to " + ae);
            }
            int nAtoms = totalEnergy.distanceMatrix().molecularSystem.size();
            header += String.format(" %12s %12s", "energy", "avgEnergy");
            values += String.format(" %12.3f %12.3f", energy, energy / nAtoms);
            for (Iterator terms = totalEnergy.energyTerms().iterator(); terms.hasNext();) {
                AbstractEnergy term = (AbstractEnergy) terms.next();
                int length = term.comment().length();
                if (length > 12) length = 12;
                header += String.format(" %12s", term.comment().substring(0, length));
                totalEnergy.offAllBut(term);
                //term.on();
                try {
                    values += String.format(" %12.3f", totalEnergy.evaluate());
                } catch (Exception ae) {
                    throw new RuntimeException("RmsGdtEnergy failed due to " + ae);
                }
            }
        }
        if (headerFlag) output.add(header);
        output.add(values);
    }


    //----------------------------------- Misc. ----------------------------------------------------

    public static void alignToX(Protein protein) {
        Atom  longestPairAtom1 = null;
        Atom  longestPairAtom2 = null;
        double longestDistance2 = -1;
        Atom atomI, atomJ;
        AtomList atoms = protein.atoms();
        int size = atoms.size();
        for (int iAtom = 0; iAtom < size; iAtom++) {
            atomI = atoms.get(iAtom);
            if (!atomI.nowhere()) {
                double iX = atomI.x();
                double iY = atomI.y();
                double iZ = atomI.z();
                for (int jAtom = iAtom+1; jAtom < size; jAtom++) {
                    atomJ = atoms.get(jAtom);
                    if (!atomJ.nowhere()) {
                        double jX = atomJ.x();
                        double jY = atomJ.y();
                        double jZ = atomJ.z();

                        double d2 = (iX - jX) * (iX - jX) + (iY - jY) * (iY - jY) + (iZ - jZ) * (iZ - jZ);
                        if (d2 > longestDistance2) {
                            longestDistance2 = d2;
                            longestPairAtom1 = atomI;
                            longestPairAtom2 = atomJ;
                        }
                    }
                }
            }
        }

        Vector3d longestVector = new Vector3d(longestPairAtom1.x()-longestPairAtom2.x(), longestPairAtom1.y()-longestPairAtom2.y(), longestPairAtom1.z()-longestPairAtom2.z());
        Vector3d xAxis = new Vector3d(1.0, 0, 0);
        double angle = xAxis.angle(longestVector);
        Vector3d rotationAxis = new Vector3d();
        rotationAxis.cross(xAxis,longestVector);
        AxisAngle4d axisAngle4d = new AxisAngle4d(rotationAxis,-angle);
        Transform3D rotationMatrix = new Transform3D();
        rotationMatrix.set(axisAngle4d);
        for (Atom atom : atoms) {
            if (!atom.nowhere()) {
                Point3d point = new Point3d(atom.x(), atom.y(), atom.z());
                rotationMatrix.transform(point);
                atom.setXYZ(point.x, point.y, point.z);
            }
        }

        if (verbose)
            System.out.println("longestAxes: "+longestPairAtom1.distanceFrom(longestPairAtom2)+"\n"+longestPairAtom1+"\n"+longestPairAtom2);

    }

    /**
     * Static version of equals. Useful when the compared objects may be null.
     * This is NOT the right place for highly sofisticate algorithms.
     */
    public static boolean equals(Object obj1, Object obj2) {
        if (obj1 == null)
            return (obj2 == null);
        else
            return (obj1.equals(obj2));
    }

    /**
     * Converts a string to an int.
     */
    public static int toInt(String s) {
        return Integer.valueOf(s.trim()).intValue();
    }

    /**
     * Converts a string to a double.
     */
    public static double toDouble(String s) {
        return Double.valueOf(s.trim()).doubleValue();
    }

    // Used for loop building

    public static int[] getSeqOfProt(Protein prot, int fromRes, int toRes) {
        int[] output = new int[toRes - fromRes + 1];
        for (int c = fromRes; c <= toRes; c++)
            if (prot.residue(c) == null)
                throw new RuntimeException("\nCurrent setting asked for type of residue " + c + "which doesn't exist in the protein\n");
            else
                output[c - fromRes] = prot.residue(c).type().ordinal();
        return output;
    }

    public static String removeBackSlashR(String sequence){
        String temp[] = sequence.split("\r");
        String out = "";
        for (String s:temp)
            out += s;
        return out;

    }


    //----------------------------------- MolecularSystem ----------------------------------------------------

    public static double radiusOfGyration(AtomList atomList) {
        ArrayList<Double> weights = new ArrayList<Double>();
        for(int i = 0; i < atomList.size(); i++) weights.add(new Double(1));
        return radiusOfGyration(atomList,weights);
    }
    public static double radiusOfGyration(AtomList atoms,ArrayList<Double> weights) {
        double sum2;
        double dx, dy, dz, d2;
        double cmx, cmy, cmz; // center of mass x, y and z
        sum2 = 0;
        cmx = cmy = cmz = 0.0;
        double nWithCoordinates = 0;
        if (atoms.size() != weights.size())
            throw new RuntimeException("This is weird.");
        for (int i = 0; i < atoms.size(); i++) {
            Atom atom = atoms.get(i);
            if (!atom.nowhere()) {
                double weight = weights.get(i).doubleValue();
                cmx += atom.x()*weight;
                cmy += atom.y()*weight;
                cmz += atom.z()*weight;
                nWithCoordinates += weight;
            }
        }
        cmx /= nWithCoordinates;
        cmy /= nWithCoordinates;
        cmz /= nWithCoordinates;
        for (Atom atom : atoms) {
            if (!atom.nowhere()) {
                dx = cmx - atom.x();
                dy = cmy - atom.y();
                dz = cmz - atom.z();
                d2 = dx * dx + dy * dy + dz * dz;
                sum2 += d2;
            }
        }
        return Math.sqrt(sum2 / nWithCoordinates);
    }

    public static void testMolecularSystemIntegrity(MolecularSystem ms, String comment) {
        int size = ms.size();
        System.out.println(comment + "   Molecular System of size = " + size);
        int normal = 0, frozen = 0, nowhere = 0, hidden = 0, clashes = 0;
        for (int i = 0; i < size; i++) {
            AtomCore ai = ms.get(i);
            if (ai.status() == AtomStatus.NORMAL) normal++;
            else if (ai.status() == AtomStatus.FROZEN) frozen++;
            else if (ai.status() == AtomStatus.NOWHERE) nowhere++;
            else if (ai.status() == AtomStatus.HIDDEN) hidden++;
            else throw new RuntimeException("MolecularSystemError:\n" + "weird AtomCore instance " + ai);
            for (int j = i + 1; j < size; j++) {
                AtomCore aj = ms.get(j);
                if (ai.atom.name.equals(aj.atom.name) &&
                        (ai.atom.residue().number() == aj.atom.residue().number()) &&
                        (ai.atom.residue().chain() == aj.atom.residue().chain()))
                    throw new RuntimeException("MolecularSystemError:\n" + "Two apparently identical getAtoms\n" + ai.atom + "\n" + aj.atom);
                if ((!ai.status().nowhere()) && (ai.x() == aj.x()) && (ai.y() == aj.y()) && (ai.z() == aj.z()))
                    throw new RuntimeException("MolecularSystemError:\n" + "Two clashing getAtoms\n" + ai.atom + "\n" + aj.atom);
            }
        }
        System.out.println("normal  = " + normal + "\n" +
                "frozen  = " + frozen + "\n" +
                "nowhere = " + nowhere + "\n" +
                "hidden  = " + hidden + "\n");
    }

    public static void hideAllBut(AtomList atomList) {
        MolecularSystem ms = atomList.molecularSystem();
        for (AtomCore atomCore : ms) {
            if ((!atomCore.atom.nowhere()) & (!atomList.contains(atomCore.atom))) atomCore.atom.hide();
        }
    }
    //----------------------------------- Optimization ----------------------------------------------------

    public static LBFGS getLBFGS(TotalEnergy energy, CommandList commands) {
        int maxSteps = commands.firstWord(MAX_STEPS).secondWordInt();
        double tolerance = commands.firstWord(TOLERANCE).secondWordDouble();
        int reportEvery = commands.firstWord(REPORT_EVERY).secondWordInt();
        try {
            return new LBFGS(energy, tolerance, maxSteps, reportEvery);
        } catch (UpdateableException ae) {
            throw new RuntimeException("getLBFGS failed due to " + ae);
        }
    }

    public static LBFGS getLBFGS(TotalEnergy energy, CommandList commands, Key key) {
        CommandList commands1 = commands.firstWordFilter(key);
        int maxSteps = commands1.secondWord(MAX_STEPS).thirdWordInt();
        double tolerance = commands1.secondWord(TOLERANCE).thirdWordDouble();
        int reportEvery = commands1.secondWord(REPORT_EVERY).thirdWordInt();
        try {
            return new LBFGS(energy, tolerance, maxSteps, reportEvery);
        } catch (UpdateableException ae) {
            throw new RuntimeException("getLBFGS failed due to " + ae);
        }
    }


    //----------------------------------- Proteins ----------------------------------------------------

    public static boolean addAtoms(Protein model,
                                   boolean allowMissingResidues,
                                   CommandList commands, PlaneParametersList parametersList, boolean fixFlag)throws UpdateableException,AlignmentException{
        //System.out.println(distanceMatrix+" 1111111111");
        boolean OK = false;
        int i = 0;
        while (!OK && i < 100000) {
            OK = addAtoms(model, allowMissingResidues, 1000, commands,parametersList, fixFlag);
            i++;
        }
        if (!OK) return false;
        return true;
    }





    public static boolean someAtomsAreNowhere(Protein model) {
        for (Atom atom : model.atoms())
            if (atom.nowhere()) return true;
        return false;
    }

    public static void resetWeirdAngles(Protein model, DistanceMatrix distanceMatrix) {
        AngleList angles = getWeirdAngles(model,  distanceMatrix);
        for (Angle angle : angles) {
            Atom atom = getLessConstraintAtom(angle);
            Utils.println("\nresetting "+atom+"\n");
            atom.resetCoordinates();
        }
    }
    public static AngleList getWeirdAngles(Protein model,DistanceMatrix distanceMatrix){
        AngleList angles = new AngleList(model.bonds(),  distanceMatrix);
        AngleList out = new AngleList();
        for (Angle angle : angles) {
            if(angle.dangerous()) out.add(angle);
        }
        return out;
    }

    private static Atom getLessConstraintAtom(Angle angle) {
        Atom out = angle.atom1();
        int n = angle.atom1().bonded().size();
        if (n > angle.atom2().bonded().size()){
            n = angle.atom2().bonded().size();
            out = angle.atom2();
        }
        if (n > angle.atom3().bonded().size()){
            n = angle.atom3().bonded().size();
            out = angle.atom3();
        }
        return out;
    }


    /*public static boolean addAtoms(Protein model,
                                   boolean allowMissingResidues,
                                   int maxNumberOfRandomCoordinatesPerResidue,
                                   CommandList commands,  String param) throws UpdateableException,AlignmentException{
        //param added
        EnergyCreator[] energyCreatorsBondedTermsOnly = {
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                new RamachandranSidechainEnergyCreator()
        };
        boolean nonDummyFound;

        if (! allowMissingResidues) {
            Residue prev;
            for (Chain chain :model.chains()){
                nonDummyFound = false;
                prev = null;
                for (Residue residue : chain) {
                    if (residue.dummy() & nonDummyFound) {
                        System.out.println("Cannot evaluate model " + model.name() + " " + residue + " is dummy");
                        return false; //there is a missing residue
                    }
                    if (!residue.dummy()) nonDummyFound = true;
                    if ((prev != null) && (residue.number() != prev.number() + 1))
                        throw new RuntimeException("Cannot add getAtoms. Missing residue(s) between " + prev + " and " + residue);
                    prev = residue;
                }
            }
        }

        model.getAtoms().freeze("In Utils.addAtoms");
        boolean OK;

        while (someAtomsAreNowhere(model)){
            DistanceMatrix dm = model.getAtoms().molecularSystem().getDistanceMatrix();
            for (Residue residue: model.residues()) {
                if (!residue.dummy()) {
					*//*for(Atom atom : residue.getAtoms().frozenAtoms()){
						if(atom.name().equals("CG")){
							atom.defrost();
						}
					}
					dm = model.getAtoms().molecularSystem().getDistanceMatrix();
					*//*

                    OK = Utils.assignRandomCoordinates(residue, commands, maxNumberOfRandomCoordinatesPerResidue);
                    dm.reset();
                    //////////////////////////tanya /////////////////////////////
                    TorsionList torsionList = createQuickAndDirtyTorsionList(residue, dm);
                    checkFrozenAtoms(torsionList);
                    dm = model.getAtoms().molecularSystem().getDistanceMatrix();
                    while (illegal(torsionList,param).size()>0) {
                        residue.getAtoms().defrostedAtoms().setNowhere();

                        OK = Utils.assignRandomCoordinates(residue, commands, maxNumberOfRandomCoordinatesPerResidue);
                        dm.reset();
                        torsionList = createQuickAndDirtyTorsionList(residue, dm);




                    if (!OK) return false;

                }
            }
            model.getAtoms().molecularSystem().createDistanceMatrix("some nowhere getAtoms got coordinates.");
            dm = model.getAtoms().molecularSystem().getDistanceMatrix();//////<----
            resetWeirdAngles(model,dm);
            for (EnergyCreator creator : energyCreatorsBondedTermsOnly) {
                creator.resetTerm();
            }
            try {
                Utils.relax(model.getAtoms(), model, energyCreatorsBondedTermsOnly, commands);
            }
            catch (EvaluationException ex) {
                dm = model.getAtoms().molecularSystem().getDistanceMatrix();
                resetWeirdAngles(model,dm);
                return false;
            }
        }
        model.getAtoms().defrost();
        return true;
    }
    }
*/

  public static boolean addAtoms(Protein model,
                                   boolean allowMissingResidues,
                                   int maxNumberOfRandomCoordinatesPerResidue,
                                   CommandList commands,  PlaneParametersList parametersList, boolean fixFlag) throws UpdateableException,AlignmentException{
        //param added
        EnergyCreator[] energyCreatorsBondedTermsOnly = {
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                new RamachandranSidechainEnergyCreator()
        };

      boolean nonDummyFound;
      if (! allowMissingResidues) {
            Residue prev;
            for (Chain chain :model.chains()){
                nonDummyFound = false;
                prev = null;
                for (Residue residue : chain) {
                    if (residue.dummy() & nonDummyFound) {
                        System.out.println("Cannot evaluate model " + model.name() + " " + residue + " is dummy");
                        return false; //there is a missing residue
                    }
                    if (!residue.dummy()) nonDummyFound = true;
                    if ((prev != null) && (residue.number() != prev.number() + 1))
                        throw new RuntimeException("Cannot add getAtoms. Missing residue(s) between " + prev + " and " + residue);
                    prev = residue;
                }
            }
        }

        model.atoms().freeze("In Utils.addAtoms");
        boolean OK;

        while (someAtomsAreNowhere(model)){
            DistanceMatrix dm ;
            dm = model.atoms().molecularSystem().getDistanceMatrix();
			Utils.assignRandomCoordinates(model.atoms(), commands, maxNumberOfRandomCoordinatesPerResidue);
            dm.reset();
                    //////////////////////////tanya /////////////////////////////
            TorsionList torsionList = createQuickAndDirtyTorsionList(model, dm);
            //checkFrozenAtoms(torsionList);
            dm = model.atoms().molecularSystem().getDistanceMatrix();
            TorsionList illegalTorsions = illegal(torsionList, parametersList);
            int counter = 0;
            while (illegalTorsions.size()>0) {
                Torsion torsion = illegalTorsions.get(0);
                Utils.println("Trying to fix "+torsion);
                Atom distantAtom = torsion.distantAtom();
                assignRandomCoordinates(distantAtom);
                dm.reset();
                if (counter++ > 100)  {
                    torsion.atoms().setNowhere();
                    torsion.atom1.residue().getAtoms().backbone().sideChains().setNowhere();
                }
                torsionList = createQuickAndDirtyTorsionList(model, dm);
                //checkFrozenAtoms(torsionList);
                illegalTorsions = illegal(torsionList, parametersList);

            }
            model.atoms().molecularSystem().createDistanceMatrix("some nowhere atoms got coordinates.");
            dm = model.atoms().molecularSystem().getDistanceMatrix();//////<----
            resetWeirdAngles(model,dm);
            for (EnergyCreator creator : energyCreatorsBondedTermsOnly) {
                creator.resetTerm();
            }
                try {
                    Utils.relax(model.atoms(), model, energyCreatorsBondedTermsOnly, commands);
                    if (fixFlag) {
                        for (Atom atom: model.atoms())
                            if (atom.energy() > 10000)
                                atom.residue().getAtoms().sideChains().setNowhere();
                    }
                }
                catch (EvaluationException ex) {
                    dm = model.atoms().molecularSystem().getDistanceMatrix();
                    dm = model.atoms().molecularSystem().getDistanceMatrix();
                    resetWeirdAngles(model,dm);
                    return false;
                }
            }
        model.atoms().defrost();
        return true;
    }

    private static TorsionList illegal(TorsionList torsionList, ParametersList parametersList) {
        TorsionList out = new TorsionList();

        for(Torsion torsion : torsionList){
            Parameters parameters = parametersList.parameters(torsion);
			/*for(Atom atom : torsion.getAtoms()){
				System.out.println(atom.frozen()+"   " + atom.name+" !!!!!!!!!!!!!!!!!");
			}

			*/
            if (parameters != null){
//                System.out.println(torsion);
//                System.out.println(parameters);
                if(createAndCheckIfElementIsLegal(torsion,parameters)==false)	out.add(torsion);
            }
        }
        Utils.println(out.size()+" illegal torsions");
        return out;

    }

    private static boolean createAndCheckIfElementIsLegal(Torsion baseElement,
                                                         Parameters parameters) {
        PlaneEnergyElement element = null;
        ////getters and setters are added
        element =  new PlaneEnergyElement(((Torsion) baseElement), parameters,0);/////<-???
        if ((element.getIsomer() == PlaneParameters.Isomer.CIS) & (element.getTorsion().cosTorsion()< 0)) {
            //throw new PlaneException(element.getTorsion(),element.getIsomer());
            Utils.println(element+" "+element.getIsomer()+"   "+element.getTorsion().cosTorsion());
            return false;
        }
        if ((element.getIsomer() == PlaneParameters.Isomer.TRANS) & (element.getTorsion().cosTorsion()>= 0)) {
            if (NQspecialCase(element) ||
                    backboneElement(element))   element.off();// This is the result of a very bad initial structure after minimization with bond and angle it will anyway get closer to the expected torsion angle.
            else {
                Utils.println(element+" "+element.getIsomer()+"   "+element.getTorsion().cosTorsion());
                return false;
            }

        }
        return true;
    }

    private static boolean backboneElement(PlaneEnergyElement element) {
        for (Atom atom : element.atoms())
            if ((!atom.isBackbone()) &&  (atom.type()!= AtomType.PCD)) return false;
        return true;
    }

    private static boolean NQspecialCase(PlaneEnergyElement element) {
        if ((element.getAtom1().type() != AtomType.NCB &&
                element.getAtom1().type() != AtomType.QCG )) return false;
        if ((element.getAtom2().type() != AtomType.NND &&
                element.getAtom2().type() != AtomType.QNE )) return false;
        if ((element.getAtom3().type() != AtomType.NCG &&
                element.getAtom3().type() != AtomType.QCD )) return false;
        if ((element.getAtom4().type() != AtomType.NOD &&
                element.getAtom4().type() != AtomType.QOE )) return false;
        return true;
    }

    private static void checkFrozenAtoms(TorsionList torsionList) {
        for(Torsion torsion : torsionList){
            int size = torsion.atoms().frozenAtoms().size();
            for(Atom atom : torsion.atoms().frozenAtoms()){
                if(size>1 && atom.name().equals("CG")){
                    atom.defrost();
                    //size--;
                }
            }

        }

    }

    private static TorsionList createQuickAndDirtyTorsionList(Protein model,
                                                              DistanceMatrix distanceMatrix) {

        AtomPairList bondList = model.bonds().filter(new GoodBonds());
        AngleList angleList = new AngleList(bondList, distanceMatrix);
        return new QuickAndDirtyTorsionList(angleList, distanceMatrix);
    }

   /* public static boolean addAtoms(Protein model,
                                   boolean allowMissingResidues,
                                   CommandList commands)throws UpdateableException,AlignmentException{
        boolean OK = false;
        int i = 0;
        while (!OK && i < 100000) {
            OK = addAtoms(model, allowMissingResidues, 1000, commands);
            i++;
        }
        if (!OK) return false;
        return true;
    }*/




/*
    public static boolean addAtoms(Protein model,
                                   boolean allowMissingResidues,
                                   int maxNumberOfRandomCoordinatesPerResidue,
                                   CommandList commands) throws UpdateableException,AlignmentException{
        EnergyCreator[] energyCreatorsBondedTermsOnly = {
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                new RamachandranSidechainEnergyCreator()
        };
        boolean nonDummyFound;

        if (! allowMissingResidues) {
            Residue prev;
            for (Chain chain :model.chains()){
                nonDummyFound = false;
                prev = null;
                for (Residue residue : chain) {
                    if (residue.dummy() & nonDummyFound) {
                        System.out.println("Cannot evaluate model " + model.name() + " " + residue + " is dummy");
                        return false; //there is a missing residue
                    }
                    if (!residue.dummy()) nonDummyFound = true;
                    if ((prev != null) && (residue.number() != prev.number() + 1))
                          throw new RuntimeException("Cannot add getAtoms. Missing residue(s) between " + prev + " and " + residue);
                    prev = residue;
                }
            }
        }

        model.getAtoms().freeze("In Utils.addAtoms");
        boolean OK;
        while (someAtomsAreNowhere(model)){
            for (Residue residue: model.residues()) {
                if (!residue.dummy()) {
                    OK = Utils.assignRandomCoordinates(residue, commands, maxNumberOfRandomCoordinatesPerResidue);
                    if (!OK) return false;

                }
            }
            model.getAtoms().molecularSystem().createDistanceMatrix("some nowhere getAtoms got coordinates.");
            DistanceMatrix dm = model.getAtoms().molecularSystem().getDistanceMatrix();
            resetWeirdAngles(model,dm);
        }
        try {
        Utils.relax(model.getAtoms(), model, energyCreatorsBondedTermsOnly, commands);
        }
        catch (EvaluationException ex) {
            DistanceMatrix dm = model.getAtoms().molecularSystem().getDistanceMatrix();
            resetWeirdAngles(model,dm);
            return false;
        }
        model.getAtoms().defrost();
        return true;
    }
  */
    public static int numberOfNonDummyResidues(Protein protein) {
        int number = 0;
        for (Iterator residues = protein.residues().iterator(); residues.hasNext();) {
            Residue residue = (Residue) residues.next();
            if (!residue.dummy()) number++;
        }
        return number;
    }

    public static int numberOfAtomsWithCoordinates(Protein protein) {
        int number = 0;
        for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (!atom.nowhere()) number++;
        }
        return number;
    }



    public static void addHydrogens(Protein protein, CommandList commands) {
        BondParametersList bondParametersList = Utils.getBondParameters(commands);
        AngleParametersList angleParametersList = Utils.getAngleParameters(commands);
        ResidueExtendedAtoms.addHydrogens(protein, bondParametersList, angleParametersList);
//	        for (Iterator iter = protein.residues().iterator(); iter.hasNext();) {
//	                Residue residue = (Residue) iter.next();
//	                if ((! residue.dummy()) & (residue instanceof ResidueExtendedAtoms))
//		                ((ResidueExtendedAtoms) residue).addHydrogens(bondParametersList,angleParametersList);
//	          }
    }

    /**
     * Creates a Protein object in a new MolecularSystem.
     */
    public static Protein getProtein(CommandList commands,
                                     Key key, ResidueCreator creator, ExceptionHandler exceptionHandler) {

        try {
            Protein out = new Protein(commands.firstWord(key).secondWord(),
                    new PdbLineATOM(),
                    creator, null);
            addHydrogens(out, commands);

            return out;
        } catch (Exception ex) {
            exceptionHandler.handle(ex);
        }

        return null;
    }

    public static Protein getProtein(CommandList commands, String fileName, ResidueCreator creator, ExceptionHandler exceptionHandler) {
        try {
            Protein out = new Protein(fileName,
                    new PdbLineATOM(),
                    creator, commands);
            addHydrogens(out, commands);
            return out;
        } catch (Exception ex) {
            exceptionHandler.handle(ex);
        }

        return null;
    }

    /**
     * Reads a protein structure from the file pointed at by the command list. The protein is created in a new MolecularSystem.
     */
    public static Protein getProtein(CommandList commands,
                                     Key key1, Key key2, ResidueCreator creator, ExceptionHandler exceptionHandler) {
        try {
            CommandList commands1 = commands.firstWordFilter(key1);
            Protein out = new Protein(commands1.secondWord(key2).thirdWord(),new PdbLineATOM(),creator, null);
            if (creator == ResidueExtendedAtoms.creator)
                addHydrogens(out, commands);
            return out;
        } catch (Exception ex) {
            exceptionHandler.handle(ex);
        }
        return null;
    }

    public static Protein getReferenceProtein(CommandList commands) {
        if (commands.keyExists(SUPERIMPOSE))
            return getProtein(commands, SUPERIMPOSE, REFERENCE, CaResidue.creator, Utils.defaultExceptionHandler);
        else return null;
    }

    public static void assignRandomCaCoordinates(Chain chain) {
        // Assign arbitrary (extended) coordinate
        Coordinates prev = new Coordinates(0.0, 0.0, 0.0);
        boolean first = true;//the dummy residue
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    residue.ca().setXYZ(prev);
                    double dx = MeshiProgram.randomNumberGenerator().nextDouble() * 3.8;
                    double tmp = Math.sqrt(3.8 * 3.8 - dx * dx);
                    double dy = MeshiProgram.randomNumberGenerator().nextDouble() * tmp;
                    double dz = Math.sqrt(3.8 * 3.8 - dx * dx - dy * dy);
                    residue.ca().addToX(dx);
                    residue.ca().addToY(dy);
                    residue.ca().addToZ(dz);
                    prev = new Coordinates(residue.ca());
                }
            }
        }
    }



    public static boolean assignRandomCoordinates(AtomList atomList, CommandList commands, int maxMissingAtoms) {
        List<Atom> nowhereAtoms = new ArrayList<Atom>();
        List<Atom> nowhereHeavyAtoms = new ArrayList<Atom>();
        List<Atom> neighbors = new ArrayList<Atom>();
        for (Atom atom : atomList) {
            if (atom.nowhere()) nowhereAtoms.add(atom);
            if (atom.nowhere() & (!atom.type().isHydrogen())) nowhereHeavyAtoms.add(atom);
        }
        if (nowhereAtoms.size() == 0) {
            println("assignRandomCaCoordinates nothing to do");
            return true;
        }
        if (nowhereHeavyAtoms.size() > maxMissingAtoms) {
            println("assignRandomCaCoordinates to missig atoms " + nowhereHeavyAtoms.size() + " maximum " + maxMissingAtoms + " do nothing");
            return false;
        }
        Utils.println("\n"+nowhereAtoms.size()+ " nowhereAtoms\n");
        for (Atom atom : nowhereAtoms) {
            Atom neighbor = null;
            for (Atom tmp : atom.bonded())
                if (!tmp.nowhere()) neighbor = tmp;
            if ((neighbor != null) && (!neighbors.contains(neighbor))){
                neighbors.add(neighbor);
            }
            else neighbors.add(null);
        }
        for (int i = 0; i < nowhereAtoms.size(); i++) {
            Atom atom = nowhereAtoms.get(i);
            Atom neighbor = neighbors.get(i);
            if (neighbor != null){
                atom.setXYZ(new Coordinates(new Coordinates(neighbor), 1.5));
                println("Assign random oordinates to " + atom+" based on "+neighbor);
            }
        }
        return true;
    }

    public static void  assignRandomCoordinates(Atom atom) {
            Atom neighbor = null;
            for (Atom tmp : atom.bonded())
                if (!tmp.nowhere()) neighbor = tmp;
            if (neighbor != null) {
                atom.setXYZ(new Coordinates(new Coordinates(neighbor), 1.5));
                println("assignRandomCaCoordinates to " + atom +" based on "+neighbor);
            }
    }

    /*
    *   as was used in Meshi.1.60
    */

    public static void assignRandomCaCoordinatesInACube(Chain chain, int boxLength, double[] bias) {
        // Assign arbitrary (extended) coordinate in a Cube of size: boxLength X boxLength X boxLength
        Coordinates prev = new Coordinates(0.0, 0.0, 0.0);
        boolean first = true;//the dummy residue
        int direction = 1;
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    residue.ca().setXYZ(prev);
                    double dx = bias[0] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double dy = bias[1] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double dz = bias[2] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double factor = direction * 3.8 * (1 / Math.sqrt(dx * dx + dy * dy + dz * dz));
                    residue.ca().addToX(dx * factor);
                    residue.ca().addToY(dy * factor);
                    residue.ca().addToZ(dz * factor);
                    prev = new Coordinates(residue.ca());
                    if (!residue.getSecondaryStructure().equals(SecondaryStructure.HELIX) && residue.position() % boxLength == 0) {
                        direction *= -1;
                        double temp = bias[0];
                        bias[0] = bias[1];
                        bias[1] = bias[2];
                        bias[2] = temp;
                    }
                }
            }
        }
    }


    public static void assignRandomCaCoordinatesInACubeRandomLength(Chain chain, int minLength, int maxLength, double[] bias) {
        // Assign arbitrary (extended) coordinate in a Cube of size: boxLength X boxLength X boxLength
        Coordinates prev = new Coordinates(0.0, 0.0, 0.0);
        boolean first = true;//the dummy residue
        int boxLength = (int) (MeshiProgram.randomNumberGenerator().nextDouble() * maxLength - minLength) + minLength;
        println("box length: " + boxLength);
        if (maxLength == minLength)
            throw new RuntimeException("the length of the box is not random");
        int direction = 1;
        boolean ss = false, changeDirection = false;
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    residue.ca().setXYZ(prev);
                    double dx = bias[0] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double dy = bias[1] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double dz = bias[2] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double factor = direction * 3.8 * (1 / Math.sqrt(dx * dx + dy * dy + dz * dz));
                    residue.ca().addToX(dx * factor);
                    residue.ca().addToY(dy * factor);
                    residue.ca().addToZ(dz * factor);
                    prev = new Coordinates(residue.ca());
                    if (residue.position() % boxLength == 0)
                        changeDirection = true;
                    if (changeDirection && residue.getSecondaryStructure().equals(SecondaryStructure.COIL)) {
                        boxLength = (int) (MeshiProgram.randomNumberGenerator().nextDouble() * maxLength - minLength) + minLength;
                        println("box length: " + boxLength);
                        direction *= -1;
                        double temp = bias[0];
                        bias[0] = bias[1];
                        bias[1] = bias[2];
                        bias[2] = temp;
                        changeDirection = false;
                    }
                }
            }
        }
    }

    public static void assignRandomCaCoordinatesInASphere(Chain chain, double sphereRadius, double[] bias) {
        // Assign arbitrary (extended) coordinate in a Cube of size: boxLength X boxLength X boxLength
        double x_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;
        double y_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;
        double z_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;

        Coordinates prev = new Coordinates(x_initial, y_initial, z_initial);
        double dis = Math.sqrt(x_initial * x_initial + y_initial * y_initial + z_initial * z_initial);
        boolean first = true;//the dummy residue
        println("dis: " + dis);
        int direction = 1;
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    residue.ca().setXYZ(prev);
                    double dx = bias[0] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double dy = bias[1] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double dz = bias[2] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    double factor = direction * 3.8 * (1 / Math.sqrt(dx * dx + dy * dy + dz * dz));
                    residue.ca().addToX(dx * factor);
                    residue.ca().addToY(dy * factor);
                    residue.ca().addToZ(dz * factor);
                    prev = new Coordinates(residue.ca());
                    dis = Math.sqrt(residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z());
                    println("dis: " + dis);
                    if (!residue.getSecondaryStructure().equals(SecondaryStructure.HELIX) && dis > sphereRadius) {
                        direction *= -1;
                        double temp = bias[0];
                        bias[0] = bias[1];
                        bias[1] = bias[2];
                        bias[2] = temp;
                    }
                }
            }
        }
    }

    public static void assignRandomCaCoordinatesInAStrongSphere(Chain chain, double sphereRadius, double[] bias) {
        // Assign arbitrary (extended) coordinate in a Cube of size: boxLength X boxLength X boxLength
        double x_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;
        double y_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;
        double z_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;

        Coordinates prev = new Coordinates(x_initial, y_initial, z_initial);
        double dis = Math.sqrt(x_initial * x_initial + y_initial * y_initial + z_initial * z_initial);
        boolean first = true;//the dummy residue
        println("good dis: " + dis);
        int direction = 1;
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    residue.ca().setXYZ(prev);
                    do {
                        double dx = bias[0] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        double dy = bias[1] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        double dz = bias[2] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        double factor = direction * 3.8 * (1 / Math.sqrt(dx * dx + dy * dy + dz * dz));
                        residue.ca().addToX(dx * factor);
                        residue.ca().addToY(dy * factor);
                        residue.ca().addToZ(dz * factor);
                        prev = new Coordinates(residue.ca());
                        dis = Math.sqrt(residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z());

                        if (!residue.getSecondaryStructure().equals(SecondaryStructure.HELIX) && dis > sphereRadius) {
                            println("bad  dis: " + dis);
                            direction *= -1;
                            double temp = bias[0];
                            bias[0] = bias[1];
                            bias[1] = bias[2];
                            bias[2] = temp;
                        } else {
                            println("good  dis: " + dis);
                        }
                    } while (!residue.getSecondaryStructure().equals(SecondaryStructure.HELIX) && dis > sphereRadius);
                }
            }
        }
    }

    public static void assignRandomCaCoordinatesInABampSphere(Chain chain, double sphereRadius) {
        double[] bias = {1, 0.25, 0.25};
        // Assign arbitrary (extended) coordinate in a Cube of size: boxLength X boxLength X boxLength
        double x_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;
        double y_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;
        double z_initial = MeshiProgram.randomNumberGenerator().nextDouble() * 2 * sphereRadius - sphereRadius;

        Coordinates prev = new Coordinates(x_initial, y_initial, z_initial);
        Coordinates tmpPrev;
        double dis = Math.sqrt(x_initial * x_initial + y_initial * y_initial + z_initial * z_initial), dis1;
        boolean first = true;//the dummy residue
        println("good dis: " + dis);
        int direction = 1;
        boolean getCloser = false;
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    if (residue.getSecondaryStructure().equals(SecondaryStructure.HELIX)) {
                        residue.ca().setXYZ(prev);
                        double dx = bias[0] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        double dy = bias[1] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        double dz = bias[2] + MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        double factor = direction * 3.8 * (1 / Math.sqrt(dx * dx + dy * dy + dz * dz));
                        residue.ca().addToX(dx * factor);
                        residue.ca().addToY(dy * factor);
                        residue.ca().addToZ(dz * factor);
                        prev = new Coordinates(residue.ca());
                    } else {
                        residue.ca().setXYZ(prev);
                        double dx = MeshiProgram.randomNumberGenerator().nextDouble() * 3.8;
                        double tmp = Math.sqrt(3.8 * 3.8 - dx * dx);
                        double dy = MeshiProgram.randomNumberGenerator().nextDouble() * tmp;
                        double dz = Math.sqrt(3.8 * 3.8 - dx * dx - dy * dy);
                        //double factor = direction * 3.8 * (1/Math.sqrt(dx*dx+dy*dy+dz*dz));
                        residue.ca().addToX(dx);//*factor);
                        residue.ca().addToY(dy);//*factor);
                        residue.ca().addToZ(dz);//*factor);
                        dis = Math.sqrt(residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z());
                        tmpPrev = new Coordinates(residue.ca());
                        for (int i = 0; i < 9; i++) {
                            residue.ca().setXYZ(prev);
                            double dx1 = MeshiProgram.randomNumberGenerator().nextDouble() * 3.8;
                            tmp = Math.sqrt(3.8 * 3.8 - dx * dx);
                            double dy1 = MeshiProgram.randomNumberGenerator().nextDouble() * tmp;
                            double dz1 = Math.sqrt(3.8 * 3.8 - dx * dx - dy * dy);
                            //factor = direction * 3.8 * (1/Math.sqrt(dx1*dx1+dy1*dy1+dz1*dz1));
                            residue.ca().addToX(dx1);//*factor);
                            residue.ca().addToY(dy1);//*factor);
                            residue.ca().addToZ(dz1);//*factor);
                            dis1 = Math.sqrt(residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z());
                            if ((getCloser && dis1 < dis) || (!getCloser && dis1 > dis)) {
                                dis = dis1;
                                tmpPrev = new Coordinates(residue.ca());
                            }
                        }
                        residue.ca().setXYZ(tmpPrev);
                        if (dis < sphereRadius / 3 || dis > sphereRadius) {
                            getCloser = !getCloser;
                            direction *= -1;
                            double temp = bias[0];
                            bias[0] = bias[1];
                            bias[1] = bias[2];
                            bias[2] = temp;
                        }
                        prev = new Coordinates(residue.ca());

                    }

                }
            }
        }
    }

    public static void assignRandomCaCoordinatesAroundASphere(Chain chain, double sphereRadius) {
        double t = 1.5;
        double x_initial = MeshiProgram.randomNumberGenerator().nextDouble() * sphereRadius;
        double tmp = Math.sqrt(sphereRadius * sphereRadius - x_initial * x_initial);
        double y_initial = MeshiProgram.randomNumberGenerator().nextDouble() * tmp;
        double z_initial = Math.sqrt(sphereRadius * sphereRadius - x_initial * x_initial - y_initial * y_initial);
        Coordinates prev = new Coordinates(x_initial, y_initial, z_initial);
        Coordinates prevPrev = new Coordinates(0, 0, 0);
        Coordinates tmpPrev;
        double dis = Math.sqrt(x_initial * x_initial + y_initial * y_initial + z_initial * z_initial), dis1;
        println("initial dis: " + dis);
        boolean first = true;//the dummy residue
        int direction = 1;
        for (Residue residue : chain) {
            if (first) first = false;
            else {
                if (!residue.dummy()) {
                    /* if (residue.getSecondaryStructure().equals(SecondaryStructure.HELIX)) {
                    residue.ca().setXYZ(prev);
                    double dx = bias[0]+MeshiProgram.randomNumberGenerator().nextDouble()*0.5;
                    double dy = bias[1]+MeshiProgram.randomNumberGenerator().nextDouble()*0.5;
                    double dz = bias[2]+MeshiProgram.randomNumberGenerator().nextDouble()*0.5;
                    double factor = direction * 3.8 * (1/Math.sqrt(dx*dx+dy*dy+dz*dz));
                    residue.ca().addToX(dx*factor);
                    residue.ca().addToY(dy*factor);
                    residue.ca().addToZ(dz*factor);
                    dis =  Math.sqrt(residue.ca().x()*residue.ca().x()+residue.ca().y()*residue.ca().y()+residue.ca().z()*residue.ca().z());
                    tmpPrev =   new Coordinates(residue.ca());
                    for (int i = 0; i < 9;i++)
                    {
                           residue.ca().setXYZ(prev);
                            double dx1 = bias[0]+MeshiProgram.randomNumberGenerator().nextDouble()*0.5;
                            double dy1 = bias[1]+MeshiProgram.randomNumberGenerator().nextDouble()*0.5;
                            double dz1 = bias[2]+MeshiProgram.randomNumberGenerator().nextDouble()*0.5;
                            factor = direction * 3.8 * (1/Math.sqrt(dx1*dx1+dy1*dy1+dz1*dz1));
                            residue.ca().addToX(dx1*factor);
                            residue.ca().addToY(dy1*factor);
                            residue.ca().addToZ(dz1*factor);
                            dis1 =  Math.sqrt(residue.ca().x()*residue.ca().x()+residue.ca().y()*residue.ca().y()+residue.ca().z()*residue.ca().z());
                            if( Math.abs(dis1-sphereRadius)<Math.abs(dis-sphereRadius)){
                                        dis = dis1;
                                        tmpPrev =   new Coordinates(residue.ca());
                            }
                    }
                    residue.ca().setXYZ(tmpPrev);
                    prev = new Coordinates(residue.ca());

                }
                else {   */
                    residue.ca().setXYZ(prev);
                    double dx, dy, dz;
                    if (Math.abs(prev.x() - prevPrev.x()) > t) {
                        dx = MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    } else {
                        dx = MeshiProgram.randomNumberGenerator().nextDouble() - 0.5;
                    }
                    if (Math.abs(prev.y() - prevPrev.y()) > t) {
                        dy = MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    } else {
                        dy = MeshiProgram.randomNumberGenerator().nextDouble() - 0.5;
                    }
                    if (Math.abs(prev.z() - prevPrev.z()) > t) {
                        dz = MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                    } else {
                        dz = MeshiProgram.randomNumberGenerator().nextDouble() - 0.5;
                    }
                    double factor = 3.8 * (1 / Math.sqrt(dx * dx + dy * dy + dz * dz));
                    residue.ca().addToX(dx * factor);
                    residue.ca().addToY(dy * factor);
                    residue.ca().addToZ(dz * factor);
                    dis = Math.sqrt(residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z());
                    println(" dis: " + dis);

                    tmpPrev = new Coordinates(residue.ca());
                    double dx1, dy1, dz1;
                    for (int i = 0; i < 9; i++) {
                        residue.ca().setXYZ(prev);
                        if (Math.abs(prev.x() - prevPrev.x()) > t) {
                            dx1 = MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                            //System.out.println(dx1);
                        } else {
                            dx1 = MeshiProgram.randomNumberGenerator().nextDouble() - 0.5;
                            // System.out.println(dx1);
                        }
                        if (Math.abs(prev.y() - prevPrev.y()) > t) {
                            dy1 = MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        } else {
                            dy1 = MeshiProgram.randomNumberGenerator().nextDouble() - 0.5;
                        }
                        if (Math.abs(prev.z() - prevPrev.z()) > t) {
                            dz1 = MeshiProgram.randomNumberGenerator().nextDouble() * 0.5;
                        } else {
                            dz1 = MeshiProgram.randomNumberGenerator().nextDouble() - 0.5;
                        }
                        factor = 3.8 * (1 / Math.sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1));
                        residue.ca().addToX(dx1 * factor);
                        residue.ca().addToY(dy1 * factor);
                        residue.ca().addToZ(dz1 * factor);
                        double dis2 = residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z();
                        println(residue.ca().x() + "," + residue.ca().y() + "," + residue.ca().z() + " dis2: " + dis2);

                        dis1 = Math.sqrt(residue.ca().x() * residue.ca().x() + residue.ca().y() * residue.ca().y() + residue.ca().z() * residue.ca().z());
                        println(" dis1 : " + dis);

                        if (Math.abs(dis1 - sphereRadius) < Math.abs(dis - sphereRadius)) {
                            dis = dis1;
                            tmpPrev = new Coordinates(residue.ca());
                        }
                    }
                    residue.ca().setXYZ(tmpPrev);
                    println("selected  dis: " + dis);
                    prevPrev = prev;
                    prev = new Coordinates(residue.ca());

                    // }

                }
            }
        }
    }

    public static Atom getCovalentNeighbor(Atom atom) {
        for (Iterator bonded = atom.bonded().iterator(); bonded.hasNext();) {
            Atom bondedAtom = (Atom) bonded.next();
            if (!bondedAtom.nowhere()) return bondedAtom;
            else if (bondedAtom.number() < atom.number())
                return getCovalentNeighbor(bondedAtom);
        }
        return getCovalentNeighborUp(atom);
    }

    /**
     * May be used only by getCovalentNeighbor(Atom atom)
     */
    private static Atom getCovalentNeighborUp(Atom atom) {
        for (Iterator bonded = atom.bonded().iterator(); bonded.hasNext();) {
            Atom bondedAtom = (Atom) bonded.next();
            if (!bondedAtom.nowhere()) return bondedAtom;
            else if (bondedAtom.number() > atom.number())
                return getCovalentNeighborUp(bondedAtom);
        }
        return null;
    }




    public static void colorByEnergy(AtomList atomList) {
        for (Atom atom : atomList) atom.resetTemperatureFactor();
        double sum = 0, sum2 = 0, avg, avg2, std, energy;
        double[] zScore;
        double min = 10000, max = -10000;


        zScore = new double[atomList.size()];
        for (Atom atom : atomList) {
            energy = atom.energy();
            sum += energy;
            sum2 += energy * energy;
        }
        avg = sum / atomList.size();
        avg2 = sum2 / atomList.size();
        std = Math.sqrt(avg2 - avg * avg);
        for (int i = 0; i < atomList.size(); i++) {
            zScore[i] = (atomList.get(i).energy() - avg) / std;
            if (zScore[i] > max) max = zScore[i];
            if (zScore[i] < min) min = zScore[i];
        }

        if ((min >= -2) & (max <= 2)) {
            colorByEnergy(atomList, min * std + avg, max * std + avg, 0, 99.9);
        } else if ((min >= -2) & (max > 2)) {
            colorByEnergy(atomList, min * std + avg, 2 * std + avg, 0, 89.9);
            colorByEnergy(atomList, 2 * std + avg, max * std + avg, 90, 99.9);
        }
        if ((min < -2) & (max <= 2)) {
            colorByEnergy(atomList, min * std + avg, -2 * std + avg, 0, 9.9);
            colorByEnergy(atomList, -2 * std + avg, max * std + avg, 10, 99.9);
        }
        if ((min < -2) & (max > 2)) {
            colorByEnergy(atomList, min * std + avg, -2 * std + avg, 0, 9.9);
            colorByEnergy(atomList, -2 * std + avg, 2 * std + avg, 10, 88.9);
            colorByEnergy(atomList, 2 * std + avg, max * std + avg, 90, 99.9);
        }
    }

    public static void colorByEnergy(AtomList atomList, double min, double max, double minColor, double maxColor) {
        double color;
        double range = max - min;
        double colorRange = maxColor - minColor;
        double energy;
        for (Atom atom : atomList) {
            energy = atom.energy();
            if ((energy >= min) & (energy <= max)) {
                color = minColor + colorRange * ((energy - min) / range);
                if (color > 99.9) {
                    throw new RuntimeException("energy = " + energy + "\n" +
                            " min & max = " + min + " " + max + "\n" +
                            "minColor & maxColor = " + minColor + " " + maxColor +
                            "range & colorRange = " + range + " " + colorRange);
                }
                try {
                    atom.setTemperatureFactor(color);
                } catch (RuntimeException ex) {
                    throw new RuntimeException(ex + "\n" + "energy = " + energy + "\n" +
                            " min & max = " + min + " " + max + "\n" +
                            "minColor & maxColor = " + minColor + " " + maxColor +
                            "range & colorRange = " + range + " " + colorRange + "\n" +
                            "color = " + color);
                }

            }
        }
    }


    public static void assignBackboneCoordinates(AtomList modelAtoms, AtomList tempAtoms) {
        Atom[] modelArray = modelAtoms.toArrayOfAtoms();
        Atom[] tempArray = tempAtoms.toArrayOfAtoms();
        Comparator comparator = new AtomComparator();
        Arrays.sort(modelArray, comparator);
        Arrays.sort(tempArray, comparator);
        int i = 0;
        Atom oxt = null;
        for (Atom modelAtom : modelArray) {
            if (!modelAtom.name().equals("OXT")) {
                while ((i < tempArray.length) &&
                        ((tempArray[i].residue().number() != modelAtom.residue().number()) |
                                (!tempArray[i].name().equals(modelAtom.name())))) {
                    i++;
                }
                if (i >= tempArray.length) throw new RuntimeException("Can't find coordinates for " + modelAtom);
                if ((!modelAtom.type().backboneCA()) ||
                        (modelAtom.distanceFrom(tempArray[i]) < 0.8))
                    modelAtom.setXYZ(new Coordinates(tempArray[i]));
                i++;
            } else oxt = modelAtom;
        }
        Atom ca = oxt.residue().ca();
        double x = ca.x() + 1;
        oxt.setXYZ(new Coordinates(x, ca.y(), ca.z()));
    }

    public static void AssignDSSP(Protein protein, CommandList commands, Key key) {
        String dsspFileName = commands.firstWord(key).secondWord();
        AssignDSSP(protein, dsspFileName);
    }

    public static void setSS(Protein model, CommandList commands) {
        if (commands.keyExists("setSS")) {
            String line = commands.firstWord("setSS").secondWord();
            String[] words = line.split(",");
            String[] numSsPair;
            int number;
            char ssChar;
            for (String word : words) {
                numSsPair = word.split("=");
                number = Integer.valueOf(numSsPair[0]);
                ssChar = numSsPair[1].charAt(0);
                for (Residue residue : model.residues()) {
                    if (residue.number() == number) {
                        residue.setSecondaryStructure(SecondaryStructure.secondaryStructure(ssChar));
                        if (Utils.verbose())
                            System.out.println("setResidue " + residue + " to " + SecondaryStructure.secondaryStructure(ssChar));
                    }
                }
            }
        } else if (Utils.verbose) System.out.println("No manual SS assignment");
    }

    public static boolean AssignDSSP(Protein protein, String dsspFileName) {
        DSSP dssp = new DSSP(dsspFileName);
        return AssignDSSP(protein, dssp);
    }

    public static boolean AssignDSSP(Protein protein, DSSP dssp) {
        if (dssp.length() == 0) return false;
        SequenceList[] OneDimViews = dssp.sequenceLists();
        MeshiSequence fullSequences[] = protein.sequences();
        for (int iChain = 0; iChain < protein.chains().size(); iChain++) {
            MeshiSequence residueSequence = OneDimViews[iChain].getResidueSequence();
            SequenceAlignment OneDimAlignment = new SequenceAlignment(OneDimViews[iChain]);

            // Hang each column isOn its first cell
            Iterator columns = OneDimAlignment.iterator();
            Iterator resIter = residueSequence.iterator();
            while (columns.hasNext()) {
                SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
                SequenceAlignmentColumn res = (SequenceAlignmentColumn) resIter.next();
                res.cell(0).addAttribute(column);
            }
            //

            SequenceAlignment alignment = SequenceAlignment.identityAlignment(fullSequences[iChain], residueSequence);
            columns = alignment.iterator();
            while (columns.hasNext()) {
                SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
                SequenceAlignmentCell cell0 = (SequenceAlignmentCell) column.cell(0);
                SequenceAlignmentCell cell1 = (SequenceAlignmentCell) column.cell(1);
                if (!cell0.gap()) {
                    Residue res = protein.chains().get(iChain).get(cell0.number);
                    if (!res.dummy()) {
                        SequenceAlignmentColumn OneDcolumn = (SequenceAlignmentColumn) cell1.getAttribute();
                        if (OneDcolumn == null) {
                            res.setSecondaryStructure('C');
                            res.setAccessibility('A');
                        } else {
                            SequenceAlignmentCell OneDcell1 = (SequenceAlignmentCell) OneDcolumn.cell(1);
                            SequenceAlignmentCell OneDcell2 = (SequenceAlignmentCell) OneDcolumn.cell(2);
                            char ss = OneDcell1.getChar();
                            char acc = OneDcell2.getChar();
                            double relativeAccesibility = Double.valueOf(OneDcell2.comment);
                            res.setSecondaryStructure(ss);
                            res.setAccessibility(acc);
                            res.setRelativeccessibility(relativeAccesibility);
                        }
                    }
                }
            }
        }
        return true;
    }

    private final static String BLANK30 = "                              ";
    public static String F(double d) {
        return String.format("\"%-10.6f\" ",d);
        //return ("\""+d+"\" "+BLANK30).substring(0,15);
    }

    public static String F(int i) {
        return ("\""+i+"\""+BLANK30).substring(0,10);
    }
    public static String F(String s) {
        return ("\""+s+"\""+BLANK30).substring(0,10);
    }
    public static boolean AssignFullDSSP(Protein protein, String dsspFileName) {
        DSSPFull dssp = new DSSPFull(dsspFileName);
        return AssignFullDSSP(protein, dssp);
    }

    public static boolean AssignFullDSSP(Protein protein, DSSP dssp) {
        if (dssp.length() == 0) return false;
        SequenceList[] OneDimViews = dssp.sequenceLists();
        MeshiSequence fullSequences[] = protein.sequences();
        for (int iChain = 0; iChain < protein.chains().size(); iChain++) {
            //Preparing the dssp sequence to be a MeshiSequence with secondary structure attributes: ss, fullSs, acc

            MeshiSequence residueSequence = OneDimViews[iChain].getResidueSequence();
            SequenceAlignment OneDimAlignment = new SequenceAlignment(OneDimViews[iChain]);

            // Hang each column isOn its first cell
            Iterator columns = OneDimAlignment.iterator();
            Iterator resIter = residueSequence.iterator();
            while (columns.hasNext()) {
                SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
                SequenceAlignmentColumn res = (SequenceAlignmentColumn) resIter.next();
                res.cell(0).addAttribute(column);
            }
            //

            SequenceAlignment alignment = SequenceAlignment.identityAlignment(fullSequences[iChain], residueSequence);
            columns = alignment.iterator();
            while (columns.hasNext()) {
                SequenceAlignmentColumn column = (SequenceAlignmentColumn) columns.next();
                SequenceAlignmentCell cell0 = (SequenceAlignmentCell) column.cell(0);
                SequenceAlignmentCell cell1 = (SequenceAlignmentCell) column.cell(1);
                if (!cell0.gap()) {
                    Residue res = protein.chains().get(iChain).get(cell0.number);
                    if (!res.dummy()) {
                        SequenceAlignmentColumn OneDcolumn = (SequenceAlignmentColumn) cell1.getAttribute();
                        if (OneDcolumn == null) {
                            res.setLocalStructure(new LocalStructure());
                            res.setSecondaryStructure('C');
                            res.setAccessibility('A');
                        } else {
                            SequenceAlignmentCell OneDcell1 = (SequenceAlignmentCell) OneDcolumn.cell(1);
                            SequenceAlignmentCell OneDcell2 = (SequenceAlignmentCell) OneDcolumn.cell(2);
                            SequenceAlignmentCell OneDcell3 = (SequenceAlignmentCell) OneDcolumn.cell(3);
                            char ss = OneDcell1.getChar();
                            char acc = OneDcell2.getChar();

                            double relativeAccesibility = Double.valueOf(OneDcell2.comment);
                            res.setSecondaryStructure(ss);
                            res.setAccessibility(acc);
                            res.setRelativeccessibility(relativeAccesibility);
                            res.setLocalStructure((LocalStructure) OneDcell3.getAttribute(MeshiAttribute.SS_DSSP));
                        }
                    }
                }
            }
        }
        return true;
    }


    public static MeshiSequence getSSPrediction(String fileName) {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String line;
        int i= 0;
        MeshiSequence ssSequence = new MeshiSequence("Sequence of SS prediction");
        try {
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("#")) continue;
                if (line.length() < 10) continue;
                ResidueSsPrediction residueSsPrediction = new ResidueSsPrediction(line);
                SequenceAlignmentCell cell = new SequenceAlignmentCell(residueSsPrediction.getResidueType().nameOneLetter().charAt(0), i);
                i++;
                ssSequence.add(cell);
                cell.addAttribute(residueSsPrediction);
            }

            reader.close();
        } catch (IOException ex){
            ex.printStackTrace();
            throw new RuntimeException(ex.getMessage());
        }
        return ssSequence;
    }

    public static AtomList getNeighbors(ResidueList residueList, Protein protein, double clashDistance) {
        List<Atom> all = new ArrayList<Atom>();
        for (Iterator residues = residueList.iterator(); residues.hasNext();) {
            Residue residue = (Residue) residues.next();
            for (Iterator atoms = residue.getAtoms().iterator(); atoms.hasNext();) {
                Atom atom = (Atom) atoms.next();
                if ((!atom.nowhere()) && (!all.contains(atom))) all.add(atom);
                List<Atom> tempList = getNeighbors(atom, protein, clashDistance);
                for (Atom candidate : tempList)
                    if ((!candidate.nowhere()) && (!all.contains(candidate))) all.add(candidate);
            }
        }
        Object[] tempArray = all.toArray();
        Arrays.sort(tempArray);

        AtomList out = new AtomList(protein.molecularSystem);
        for (Object atom : tempArray) {
            if (weirdAtom((Atom) atom))
                throw new RuntimeException("This is weird " + atom + " " + ((Atom) atom).nowhere());
            out.add((Atom) atom);
        }
        return out;
    }

    public static boolean weirdAtom(Atom atom) {
        if (atom.nowhere()) return true;
        if (atom.x() < -1000) return true;
        if (atom.y() < -1000) return true;
        if (atom.z() < -1000) return true;
        if (atom.x() > 1000) return true;
        if (atom.y() > 1000) return true;
        if (atom.z() > 1000) return true;
        return false;
    }

    public static List<Atom> getNeighbors(Atom atom, Protein protein, double clashDistance) {
        ArrayList<Atom> out = new ArrayList<Atom>();
        for (Iterator atoms = protein.atoms().iterator(); atoms.hasNext();) {
            Atom other = (Atom) atoms.next();
            if ((!atom.nowhere()) && (atom != other) && (!other.nowhere()))
                if (atom.distanceFrom(other) < clashDistance * 3)
                    out.add(other);
        }
        return out;
    }

    public static int getClashes(Atom atom, List<Atom> neighbors, double clashDistance) {
        int out = 0;
        for (Atom neighbor : neighbors) {
            if (atom.distanceFrom(neighbor) < clashDistance) out++;
            if (atom.distanceFrom(neighbor) < 2) out += 20;
        }
        return out;
    }


    //----------------------------------- Sequences ----------------------------------------------------

    public static ResidueSequence getResidueSequence(CommandList commands, Key key) {
        String fileName = commands.firstWord(key).secondWord();
        FastaList fastaList = new FastaList(fileName);
        return new ResidueSequence(fastaList);
    }


    //----------------------------------- Lists ----------------------------------------------------

    public static void print(Iterable list) {
        for (Object obj : list)
            System.out.println(obj);
    }

    /**
     * Print the list elements in raws.
     * The number of elements in a raw is determined by the parameter
     */
    public static void print(ArrayList list, int rawLength, String format) {
        if (rawLength < 1) throw new RuntimeException("Raw Length cannot be shorter than one.");
        int i = 0;
        for (Object obj : list) {
            System.out.printf(format, obj);
            if (i % rawLength == 0) System.out.println();
            i++;
        }
        System.out.println();
    }

    public static ArrayList filter(ArrayList list, Filter filter, ArrayList out) {
        if (out.size() != 0) throw new RuntimeException("out needs to be empty");
        for (Object obj : list)
            if (filter.accept(obj)) out.add(obj);
        return out;
    }

    public static DistanceLists filter(DistanceLists list, Filter filter, DistanceLists out) {
        if (out.size() != 0) throw new RuntimeException("out needs to be empty");
        for (DistanceList distanceRow : list) {
            DistanceList newRow = new DistanceList(distanceRow.atomOne);
            out.add(newRow);
            for (Distance dis : distanceRow)
                if (filter.accept(dis)) newRow.add(dis);
        }
        return out;
    }

    //-------------------------------------------------------------------------------------------------------------------------------------------

    public static void testNonBondedList(DistanceLists nonBondedList) {
        FreeDistance freeDistance;
        for (DistanceList distanceRow : nonBondedList) {
            for (Distance distance : distanceRow) {
                if (distance.mode() == DistanceMode.INFINITE)
                    throw new RuntimeException("Infinite distance in NonBondedList" + distance);
                freeDistance = new FreeDistance(distance.atom1(), distance.atom2());
                if (freeDistance.distance() != distance.distance())
                    throw new RuntimeException("Unupdated distance in NonBondedList:" + distance + "\nReal distance is " + freeDistance);
                if (distance.distance() >= 5.5)
                    throw new RuntimeException("distance > 5.5 "+distance);
            }
        }
    }
    //------------------------------------------------------------------------------------------------------------------------------------------   

    public static AtomList duplicateInAnewMolecularSystem(AtomList atomList, MolecularSystem newMolecularSystem) {
        AtomList out = new AtomList(newMolecularSystem);

        for (Atom atom : atomList) {
            //    public Atom(String name, Residue residue, AtomType type, Coordinates coordinates, double occupency, Double temperatureFactor,MolecularSystem molecularSystem) {
            Residue residue = atom.residue();
            Residue newResidue = new Residue(new ResidueIdentifier(residue.chain(), residue.getChainNumber(),
                                                                   residue.number()),
                                             residue.name);
            Coordinates coordinates = new Coordinates(atom.x(), atom.y(), atom.z());
            out.add(new Atom(atom.name, newResidue, atom.type(), coordinates, atom.occupancy(), atom.temperatureFactor(), newMolecularSystem));
        }
        return out;
    }

    public static void moveHydrogensCloseToHeavy(AtomList atoms, double radius) {
        for (Iterator atomsIter = atoms.iterator(); atomsIter.hasNext();) {
            Atom atom = (Atom) atomsIter.next();
            if (atom.type().isHydrogen()) {
                Atom heavy = atom.bonded().get(0);
                atom.randomize(radius, heavy);
            }
        }
    }

    public static BondParametersList getBondParameters(CommandList commands) {
        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        return new BondParametersList(parametersDirectory +
                "/" + BOND_PARAMETERS);
    }

    public static AngleParametersList getAngleParameters(CommandList commands) {
        Command command = commands.firstWord(PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        return new AngleParametersList(parametersDirectory +
                "/" + ANGLE_PARAMETERS);
    }


    private static class DefaultExceptionHandler implements ExceptionHandler {
        public void handle(Exception ex) {
            System.out.println("An Exception has occured and handled by the  DefaultExceptionHandler.\n" + ex);
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    private static class DoNothingExceptionHandler implements ExceptionHandler {
        public void handle(Exception ex) {
            System.out.println("An Exception " + ex + " has occured.\n" + "The doNothing ExceptionHandler did nothing.");
        }
    }

    //---------------------------------------------------------------------- STRICT -------------------------------------------------
    private static boolean strict = true;
    public static void setStrictOn() {
        strict = true;
        Utils.println("Strict is ON");
    }
    public static void setStrictOff() {
        strict = false;
        Utils.println("Strict is OFF");
    }
    public static boolean isStrict() {
        return strict;
    }

    //---------------------------------------------------------------------- VERBOSE -------------------------------------------------
    private static boolean verbose = true;

    public static boolean verbose() {
        return verbose;
    }

    public static void verboseOn() {
        Utils.println("Verbose mode ON");
        verbose = true;
    }

    public static void verboseOff() {
        Utils.println("Verbose mode OFF");
        verbose = false;
    }

    public static void println(String s) {
        if (verbose) System.out.println(s);
    }
     public static void throwException(Object obj,Exception ex, String comment) {
        System.out.println("******************** ERROR in "+obj+" *************************");
        System.out.println(ex.getMessage());
        System.out.println(comment);
        ex.printStackTrace();
        throw new RuntimeException("Quiting");
    }

    public static void println() {
        if (verbose) System.out.println();
    }

    public static void print(String s) {
        if (verbose) System.out.print(s);
    }


    public static void printDebug(Object obj,String s) {
        System.out.println("Debug: "+obj+" ("+s+")");
    }
    public static void printDebug(Object obj,Protein protein, String fileName) {
        Utils.printDebug(obj, "******** debug ********** Printing debug file " + fileName);
        try {
            MeshiWriter writer = new MeshiWriter("debug." + fileName + ".pdb");
            protein.atoms().print(writer);
            writer.close();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        Utils.printDebug(obj, "******** closing debug **********");
    }
    //-----------------------------------------------------------------------------------

    public static ArrayList<RamachandranEnergyElement> ramachEnergyPerResidue(RamachandranEnergy ramach) {
        ArrayList<RamachandranEnergyElement> out = new ArrayList<RamachandranEnergyElement>();
        for (Object o : ramach.elementsList()) {
            RamachandranEnergyElement element = (RamachandranEnergyElement) o;
            double e = element.evaluate() / element.weight();
            //if (e > 5) Utils.println("Ramachandran element "+element+" ; energy = "+e);
            if (e > RamachandranEnergyElement.THRESHOLD) out.add(element);
        }
        return out;
    }

    //----------------------- miscellaneous ------------------------
    public static boolean isAbnormal(double x) {
        if (x < 0) return false;
        if (x == 0) return false;
        if (x > 0) return false;
        return true;
    }

    public static AtomPairList detectContacts(AtomList atoms, Filter filter1, Filter filter2,double threshold,AtomPairList outList) {
        if (outList.size() > 0)
            throw new RuntimeException("outList must be empty");
        int size  = atoms.size();
        for (int iAtom = 0; iAtom < size; iAtom++) {
            Atom atomI = atoms.get(iAtom);
            if (filter1.accept(atomI)) {
                for (int jAtom = iAtom+1; jAtom < size; jAtom++){
                    Atom atomJ = atoms.get(jAtom);
                    if (filter2.accept(atomJ)){
                        double d = atomI.distanceFrom(atomJ);
                        if (d < threshold) outList.add(new AtomPair(atomI,atomJ));
                    }

                }
            }
        }
        return outList;
    }

    public static AtomPairList detectContacts(AtomList atoms, Filter filter, double threshold) {
      return detectContacts(atoms,filter,filter,threshold,new AtomPairList());
    }

    public static String pdb2Fasta(File in) throws IOException{
        AtomList atoms = new AtomList(in.getAbsolutePath());
        new MolecularSystem();
        Protein protein = new Protein(atoms, ResidueExtendedAtomsCreator.creator);
        String fastaName = protein.name()+".fasta";
        MeshiWriter writer = new MeshiWriter(fastaName);
        for (Chain chain : protein.chains()) {
            System.out.println(">" + protein.name() + " " + chain.name() + " ; " + chain.toString());
            writer.println(">" + protein.name() + " " + chain.name() + " ; " + chain.toString());
        }
        writer.close();
        return fastaName;
    }
}
