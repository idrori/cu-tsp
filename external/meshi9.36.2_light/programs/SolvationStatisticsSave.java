package programs;

import meshi.energy.contacts.ContactsAndSASA;
import meshi.energy.contacts.ContactsCreator;
import meshi.energy.contacts.ContactsInfo;
import meshi.energy.solvation.SolvationCreator;
import meshi.energy.solvation.SolvationEnergy;
import meshi.energy.solvation.SolvationInfoElement;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomCore;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.AtomType;
import meshi.parameters.SecondaryStructure;
import meshi.sequences.AlignmentException;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.mathTools.Spline2D;

import java.io.File;
import java.io.IOException;

/**
 */
public class SolvationStatisticsSave {
    //public static final double BB_POLAR_WEIGHT = 4.2001, BB_NON_POLAR_WEIGHT = 0.3312, SC_POLAR_WEIGHT =  3.3619, SC_NON_POLAR_WEIGHT = 0.3939;
    public static final double BB_POLAR_WEIGHT = 1, BB_NON_POLAR_WEIGHT = 1, SC_POLAR_WEIGHT =  1, SC_NON_POLAR_WEIGHT =1;
    private double maxCNC = -1;
    private double maxHBC = -1;
    private int[][] cncHistograms = new int[AtomType.values().length][80];
    private int[][][] hbcHistograms = new int[AtomType.values().length][80][80];
    private String nativeFileName = null;
    private double[] gdt = {0, 0, 0, 0, 0};
    private File outputDir;
    private String targetName = null;
    private MeshiWriter targetSolvationData = null, modelNames = null;
    private int modelNumber = -1;
    private String oldNativeFileName = null;
    private Protein nativeStructure = null;
    private boolean nativeFlag;
    private File targetDir = null;
    private File modelsDir = null;
    private File modelsOutputDir = null;
    public static void main(String[] args) throws IOException, AlignmentException, UpdateableException {
        MeshiProgram.initRandom(0);
        SolvationStatisticsSave solvationStatistics = new SolvationStatisticsSave();
        solvationStatistics.collectStatistics(args[0], args[1], args[2]);
    }

   public SolvationStatisticsSave() {}


    public void initHistograms() {
        for (int iType = 0; iType < AtomType.values().length; iType++)
            for (int jCNC = 0; jCNC < 30; jCNC++) {
                cncHistograms[iType][jCNC] = 0;
                for (int kHBC = 0; kHBC < 8; kHBC++)
                    hbcHistograms[iType][jCNC][kHBC] = 0;
            }
    }

//    if (!nativeFlag) {
//        //Utils.printDebug(this,outputDir.getPath());
//        modelsDir = new File(outputDir.getPath()+"/models");
//        modelsDir.mkdir();
//        targetSolvationData = new MeshiWriter(outputDir.getPath()+"/energies.dat");
//        modelNames = new MeshiWriter(outputDir.getPath()+"/names.txt");
//    }

    public void initTarget() throws IOException{
        targetDir = new File(outputDir.getPath()+"/"+targetName);
        if (targetDir.exists()) throw new RuntimeException("This is weird. "+targetDir);
        targetDir.mkdir();
        modelsDir = new File(targetDir.getPath()+"/models");
        modelsDir.mkdir();

        if (modelNames != null) modelNames.close();
        String targetNamesFileName =  targetDir.getPath()+"\\"+targetName+"_names.txt";
        System.out.println("targetNamesFileName = "+targetNamesFileName);
        modelNames = new MeshiWriter(targetNamesFileName);
        if (targetSolvationData != null) targetSolvationData.close();
        targetSolvationData = new MeshiWriter(targetDir.getPath()+"/"+targetName+"_energies.dat");
        modelNumber = 0;
        nativeStructure = new Protein(new AtomList(nativeFileName), new ResidueExtendedAtomsCreator());
    }

    public void analyzeModel(String fileName, SolvationCreator solvationCreator, ContactsCreator contactsCreator, CommandList commands) throws AlignmentException, UpdateableException, IOException{
        Protein model;
        new MolecularSystem();
        model = new Protein(new AtomList(fileName), new ResidueExtendedAtomsCreator());
        for (Atom atom : model.atoms()) {
            if((!atom.nowhere()) && (atom.x() < -999))
                  throw new RuntimeException("This is weird.\n"+atom+"  "+atom.x()+" "+(atom.x() < -999));
        }
        Utils.addHydrogens(model,commands);
        if (nativeStructure != null)
            gdt = Rms.gdt(nativeStructure, model);
        removeOutliers(model);

        String dsspFileName;
        if (!nativeFlag)
            dsspFileName = "dssp/" + targetName+"/"+model.name() + ".dssp";
        else dsspFileName = "dssp/" + model.name() + ".dssp";
        if (!Utils.AssignDSSP(model, dsspFileName)) return;
        for (int iResidue = 0; iResidue < model.residues().size(); iResidue++) {
            Residue residue = model.residues().get(iResidue);
            if (residue.getSecondaryStructure() == SecondaryStructure.UNK) {
                residue.setSecondaryStructure(SecondaryStructure.COIL);
            }
        }

        SolvationEnergy solvationEnergy;
        ContactsAndSASA contactsAndSASA = (ContactsAndSASA) contactsCreator.createEnergyTerm(model,null,commands);
        solvationEnergy = (SolvationEnergy) solvationCreator.createEnergyTerm(model, model.atoms().molecularSystem.getDistanceMatrix(), commands);

        solvationEnergy.update(1);
        SolvationInfoElement solvationInfoElement = (SolvationInfoElement) solvationEnergy.evaluate();
        ContactsInfo         contactsInfo         = (ContactsInfo)          contactsAndSASA.evaluate();

        double[] cnc = solvationEnergy.cnc();
        double[] hbc = solvationEnergy.hbc();
        double[] energies = solvationEnergy.atomEnergies;

        MeshiWriter modelData = null;
        modelData = new MeshiWriter(modelsDir.getPath()+"/"+model.name()+"_energies.dat");
        MolecularSystem molecularSystem = model.atoms().molecularSystem;
        for (int iAtom = 0; iAtom < molecularSystem.size(); iAtom++) {
            AtomCore atomCore = molecularSystem.get(iAtom);
           //if (atomCore.atom.residue().getSecondaryStructure().equals(SecondaryStructure.COIL)) continue;
           modelData.println(atomCore.number+"\t"+atomCore.type()+"\t"+atomCore.atom.name()+"\t"+atomCore.atom.residueNumber()+"\t"+
                               cnc[iAtom]+"\t"+hbc[iAtom]+"\t"+energies[iAtom]+"\t"+atomCore.atom.nowhere());
            AtomType type = atomCore.type();
            int ordinal = type.ordinal();
                    int cncIndex = (int) Math.floor(cnc[iAtom] - 0.5) + 1;
                    cncHistograms[ordinal][cncIndex]++;
                    if (cnc[iAtom] > maxCNC) maxCNC = cnc[iAtom];
                    int hbcIndex = (int) Math.floor(hbc[iAtom] - 0.5) + 1;
                    hbcHistograms[ordinal][cncIndex][hbcIndex]++;
                    if (hbc[iAtom] > maxHBC) maxHBC = hbc[iAtom];
        }


            modelNames.println(fileName + "  " + modelNumber);
            double[] atomTypeEnergies = new double[AtomType.values().length];
            for (int i = 1; i< AtomType.values().length; i++)
                atomTypeEnergies[i] = 0;
            for (int i = 0; i < model.atoms().size(); i++) {
                Atom atom = model.atoms().get(i);
                atomTypeEnergies[atom.type().ordinal()] += energies[atom.number()];
            }

            targetSolvationData.print(modelNumber + "\t" + model.chain().getNonDummyResidues().size() + "\t" + gdt[0] + "   " + solvationInfoElement.getValue() + "\t" + "\t");
            targetSolvationData.print(solvationInfoElement.bbPolarInfo.getValue() + "\t" + solvationInfoElement.bbCarbonInfo.getValue() + "\t" + solvationInfoElement.scPolarInfo.getValue() + "\t" + solvationInfoElement.scCarbonInfo.getValue()+"\t"+solvationInfoElement.hbInfo.getValue()+"\t"+solvationInfoElement.stdInfo.getValue()+"\t"+solvationInfoElement.entropyInfo.getValue()+"\t");
            targetSolvationData.print(solvationInfoElement.bbPolarN.getValue() + "\t" + solvationInfoElement.bbCarbonN.getValue() + "\t" + solvationInfoElement.scPolarN.getValue() + "\t" + solvationInfoElement.scCarbonN.getValue());
            targetSolvationData.print("\t"+contactsInfo.contacts12.getValue()+"\t"+contactsInfo.contacts14.getValue()+"\t"+contactsInfo.contacts14CoreRatio.getValue()+"\t"+contactsInfo.sasaRatio.getValue());
        for (int i = 0; i< AtomType.values().length; i++) {
                targetSolvationData.print("\t" + atomTypeEnergies[i]);
            }
            targetSolvationData.println();
            modelNumber++;
        modelData.close();
    }
    public void collectStatistics(String fileNames, String outputDirectoryName, String nativeFlagString) throws IOException, AlignmentException, UpdateableException{
        String fileName;
        Runtime runtime = Runtime.getRuntime();
        SolvationCreator solvationCreator = new SolvationCreator(1.0, SC_POLAR_WEIGHT, SC_NON_POLAR_WEIGHT, BB_POLAR_WEIGHT, BB_NON_POLAR_WEIGHT, 0.0001);
        ContactsCreator  contactsCreator  = new ContactsCreator();
        nativeFlag = nativeFlagString.equals("native");
        if ((!nativeFlag) & (!nativeFlagString.equals("decoys")))
            throw new RuntimeException("nativeFlagString my either be \"native\" or \"decoys\"");

        MeshiLineReader fileNamesReader = new MeshiLineReader(fileNames);
        MeshiWriter proteinData;
        outputDir = new File(outputDirectoryName);
        if (outputDir.exists())
            throw new RuntimeException("Please remove old output dir");
        outputDir.mkdir();
        if (nativeFlag) {
            initHistograms();
            modelsOutputDir = new File(outputDir.getPath()+"/models");
            modelsOutputDir.mkdir();
        }

        CommandList commands = new CommandList("commands.MCM");
        boolean first = true;
        while ((fileName = fileNamesReader.readLine()) != null) {
            System.out.println("Analyzing "+fileName);
            if (!nativeFlag) {
                nativeFileName = fileName.substring(0, 36) + "N.pdb";
                targetName = fileName.substring(27,35);
                System.out.println("TargetName = "+targetName+" Native file name "+nativeFileName+" oldNativeFileName "+oldNativeFileName);
                if (!nativeFileName.equals(oldNativeFileName)) {//that is a new target
                    oldNativeFileName = nativeFileName;
                    initTarget();
                }
            } else {
                if (first) {
                    modelsDir = new File(outputDir.getPath() + "/models");
                    modelNames = new MeshiWriter(modelsDir.getAbsolutePath() + "\\modelsNames.txt");
                    targetSolvationData = new MeshiWriter(modelsDir.getAbsolutePath() + "\\modelsEnergies.dat");
                    first = false;
                }
            }
            analyzeModel(fileName, solvationCreator, contactsCreator, commands );
            runtime.gc();
        }
        if (nativeFlag) {
            writeXML(outputDir.getPath()+"/CNCdata.xml", "CNC", cncHistograms, null, "CNC = Carbon Neighbor Count");
            writeXML(outputDir.getPath()+"/HBCdat.xml" , "HBC", null, hbcHistograms, "HBC = Hydrogen Bond Count");
        }
        System.out.println("maxCNC = " + maxCNC);
        System.out.println("maxHBC = " + maxHBC);
        targetSolvationData.close();
        modelNames.close();
    }
    public static void writeXML(String fileName, String label, int[][] histograms1D, int[][][] histograms2D, String comment) throws IOException{
        MeshiWriter writer = new MeshiWriter(fileName+label+".xml");

        writer.println("<?xml version=\"1.0\" encoding=\"utf-8\"?>");
        //writer.println("<comment>"+comment+"</comment>");
        writer.println();
        writer.println("<"+label+"parameters>");
        for (AtomType type : AtomType.values()) {
            if ((type == AtomType.XXX )|| (type == AtomType.HE)) continue;
            int ordinal = type.ordinal();
            String group, backbone;
            if (!type.isCarbon())
                if (type == AtomType.PN)
                    group = "nonPolar";
                else
                    group = "polar";
            else group = "nonPolar";
            if (type.backbone())
                backbone = "backbone";
            else backbone = "sideChain";
            writer.print("\t<"+label+" name=\""+type+"\" group=\""+group+"\" backbone=\""+backbone+"\">");
            if (histograms1D != null) {
                int last = -1;
                for (int iBinCNC = 0; iBinCNC < histograms1D[0].length; iBinCNC++)
                    if (histograms1D[ordinal][iBinCNC] > 0) last = iBinCNC;
                if (last != 0)
                    for (int iBinCNC = 0; iBinCNC <= last; iBinCNC++) {
                        writer.print(histograms1D[ordinal][iBinCNC] + " ");
                    }
            }
            else {
                int lastCNC = -1, lastHBC = -1;
                for (int iBinCNC = 0; iBinCNC < histograms2D[0].length; iBinCNC++)
                    for (int jBinHBC = 0; jBinHBC < histograms2D[0][0].length; jBinHBC++)
                    if (histograms2D[ordinal][iBinCNC][jBinHBC] > 0) {
                        if (lastCNC < iBinCNC) lastCNC = iBinCNC;
                        if (lastHBC < jBinHBC) lastHBC = jBinHBC;
                    }
                if (lastCNC != 0)
                    for (int iBinCNC = 0; iBinCNC <= lastCNC; iBinCNC++) {
                        writer.print("\n\t"+"<CNCbin>"+iBinCNC);
                        for (int jBinHBC = 0; jBinHBC <= lastHBC; jBinHBC++)
                            writer.print("<HBCbin>  "+jBinHBC+" "+histograms2D[ordinal][iBinCNC][jBinHBC]+"</HBCbin>\t");
                        writer.print("</CNCbin>");
                    }
            }
            writer.println("</"+label+">");
        }
        writer.println("</" + label + "parameters>");
        writer.close();
    }


    public static void removeOutliers(Protein protein) {
        for (Atom atom:protein.atoms())  {
            if ((!atom.nowhere()) && (atom.x() < -999))
                atom.setNowhere();
        }
    }

    public static void testProtein(Protein protein,SolvationEnergy solvationEnergy) throws IOException{
        MeshiWriter writer = new MeshiWriter("proteinTest.dat");
        solvationEnergy.evaluate();
        double[] cnc = solvationEnergy.cnc();
        double[] hbc = solvationEnergy.hbc();
        double[] energies = solvationEnergy.atomEnergies;
        MolecularSystem molecularSystem = protein.atoms().molecularSystem;
        testAtom(writer, molecularSystem, cnc, hbc, solvationEnergy, 1);
        for (int iAtom = 0; iAtom < molecularSystem.size(); iAtom++) {
            writer.println(molecularSystem.get(iAtom).type()+", "+molecularSystem.get(iAtom).atom + "  " + cnc[iAtom] + "   " + hbc[iAtom] + "   " + energies[iAtom]);
        }
        writer.close();
    }

    static void testAtom(MeshiWriter writer, MolecularSystem molecularSystem, double[] CNC, double[] HBC, SolvationEnergy solvationEnergy, int atomNumber){
        double cnc = CNC[atomNumber];
        double hbc = HBC[atomNumber];
        AtomCore atom = molecularSystem.get(atomNumber);
        Spline2D spline2D = solvationEnergy.splines()[atomNumber];
        writer.println("Testing "+cnc+" , "+hbc+" , "+spline2D+" , "+atom);
        spline2D.calc(cnc,hbc);
        writer.println("spline2D.s = "+spline2D.s);
        spline2D.test(writer, cnc, hbc);
    }

}
