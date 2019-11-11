package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EnergyInfoElement;
import meshi.energy.TotalEnergy;
import meshi.energy.solvation.*;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;
import meshi.parameters.SecondaryStructure;
import meshi.util.*;
import meshi.util.file.MeshiWriter;
import meshi.util.info.*;
import meshi.util.statistics.*;


import java.io.IOException;
import java.util.ArrayList;

/**
 */
public class SolvationProject {

    public static void main(String[] args) {
        ProteinAnalyzer[] analyzers;
        Utils.verboseOff();
        try {
            MeshiProgram.initRandom(0);
            CommandList commands = new CommandList(args[0]);
            SolvationAnalyzer solvationAnalyzer = new SolvationAnalyzer(commands);
            if (args[3].equals("NATIVE")) {
                analyzers = new ProteinAnalyzer[1];
                analyzers[0] = solvationAnalyzer;
            }
            else {
                analyzers = new ProteinAnalyzer[3];
                analyzers[0] = new GdtRmsAnalyzer(); analyzers[1] = solvationAnalyzer; analyzers[2] = new GoapStatistics(commands);
            }
            GetStataistics.getStatistics(args, analyzers);
            SolvationAnalyzer.writeXML(args[2] + "\\SolvationStatistics.", "All",   solvationAnalyzer.hbcHistogramsAll,   solvationAnalyzer.nNowhereAll);
            SolvationAnalyzer.writeXML(args[2] + "\\SolvationStatistics.", "Helix", solvationAnalyzer.hbcHistogramsHelix, solvationAnalyzer.nNowhereHelix);
            SolvationAnalyzer.writeXML(args[2] + "\\SolvationStatistics.", "Sheet", solvationAnalyzer.hbcHistogramsSheet, solvationAnalyzer.nNowhereSheet);
            SolvationAnalyzer.writeXML(args[2] + "\\SolvationStatistics.", "Coil",  solvationAnalyzer.hbcHistogramsCoil,  solvationAnalyzer.nNowhereCoil);

        } catch (Exception ex) {
            ex.printStackTrace();
            throw new RuntimeException(ex);
        }
    }

    static class SolvationAnalyzer implements ProteinAnalyzer {
        CommandList commands;
        int nTypes = AtomType.values().length;
        int[][][] hbcHistogramsAll = new int[nTypes][80][80];
        int[][][] hbcHistogramsCoil = new int[nTypes][80][80];
        int[][][] hbcHistogramsHelix = new int[nTypes][80][80];
        int[][][] hbcHistogramsSheet = new int[nTypes][80][80];
        int[] nNowhereAll = new int[nTypes];
        int[] nNowhereHelix = new int[nTypes];
        int[] nNowhereSheet = new int[nTypes];
        int[] nNowhereCoil = new int[nTypes];
        EnergyCreator[] creator = {new SolvationCreator()};

        public SolvationAnalyzer(CommandList commands) {
            this.commands = commands;
            EnergyCreator[] creator = {new SolvationCreator()};
            initHistogram(hbcHistogramsAll);
            initHistogram(hbcHistogramsCoil);
            initHistogram(hbcHistogramsHelix);
            initHistogram(hbcHistogramsSheet);
            for (int i = 0; i < nTypes; i++) {
                nNowhereAll[i] = nNowhereHelix[i] = nNowhereSheet[i] = nNowhereCoil[i] = 0;
            }
        }

        public static void initHistogram(int[][][] hbcHistograms) {
            for (int iType = 0; iType < hbcHistograms.length; iType++)
                for (int jCNC = 0; jCNC < hbcHistograms[0].length; jCNC++) {
                    for (int kHBC = 0; kHBC < hbcHistograms[0][0].length; kHBC++)
                        hbcHistograms[iType][jCNC][kHBC] = 0;
                }
        }

        public ProteinInfo analyze(Protein protein) {
            try {
                ArrayList<MeshiInfo> infoList = new ArrayList();
                ProteinInfo info;
                EnergyInfoElement infoElement;
                SolvationEnergy solvationEnergy;
                TotalEnergy energy = new TotalEnergy(protein, creator, commands, "");
                energy.update();
                energy.evaluate();
                solvationEnergy = (SolvationEnergy) energy.getEnergyTerm(new SolvationEnergy());
                infoElement = solvationEnergy.evaluate();
                infoList.add(infoElement);
                info = new ProteinInfo(protein.metaData(), infoList, protein);
                double[] cnc = solvationEnergy.cnc();
                double[] hbc = solvationEnergy.hbc();


                for (Atom atom : protein.atoms()) {
                    if (atom.isHydrogen()) continue;
                    AtomType type = atom.type();
                    int ordinal = type.ordinal();
                    int cncIndex = (int) Math.floor(cnc[atom.number()] - 0.5) + 1;
                    int hbcIndex = (int) Math.floor(2*hbc[atom.number()] - 0.5) + 1;            //19.1.16
                    if (atom.nowhere()) nNowhereAll[ordinal]++;
                    else hbcHistogramsAll[ordinal][cncIndex][hbcIndex]++;
                    SecondaryStructure secondaryStructure = atom.residue().getSecondaryStructure();
                    switch (secondaryStructure) {
                        case HELIX:
                            if (atom.nowhere()) nNowhereHelix[ordinal]++;
                            else hbcHistogramsHelix[ordinal][cncIndex][hbcIndex]++;
                            break;
                        case SHEET:
                            if (atom.nowhere()) nNowhereSheet[ordinal]++;
                            else hbcHistogramsSheet[ordinal][cncIndex][hbcIndex]++;
                            break;
                        default:
                            if (atom.nowhere()) nNowhereCoil[ordinal]++;
                            else hbcHistogramsCoil[ordinal][cncIndex][hbcIndex]++;
                            break;
                    }
                }
                return info;
            } catch (Exception ex) {
                ex.printStackTrace();
                throw new RuntimeException(ex);
            }
        }

        public static void writeXML(String fileName, String label, int[][][] histograms2D, int[] nNowhereAtoms) throws IOException {
            MeshiWriter writer = new MeshiWriter(fileName + label + ".xml");

            writer.println("<?xml version=\"1.0\" encoding=\"utf-8\"?>");
            writer.println();
            writer.println("<" + label + "_statistics>");
            for (AtomType type : AtomType.values()) {
                if ((type == AtomType.XXX) || (type == AtomType.HE)) continue;
                int ordinal = type.ordinal();
                int lastCNC = -1, lastHBC = -1;
                for (int iBinCNC = 0; iBinCNC < histograms2D[0].length; iBinCNC++)
                    for (int jBinHBC = 0; jBinHBC < histograms2D[0][0].length; jBinHBC++)
                        if (histograms2D[ordinal][iBinCNC][jBinHBC] > 0) {
                            if (lastCNC < iBinCNC) lastCNC = iBinCNC;
                            if (lastHBC < jBinHBC) lastHBC = jBinHBC;
                        }

                String group, backbone;
                if (type.backbone())
                    backbone = "backbone";
                else backbone = "sideChain";

                if (lastCNC == 0) {
                    group = "other";
                    writer.print("\t<" + label + " name=\"" + type + "\" group=\"" + group + "\" backbone=\"" + backbone + "\" nNowhere=\"0\"/>\n");
                } else {
                    if (!type.isCarbon())
                        if (type == AtomType.PN)
                            group = "nonPolar";
                        else
                            group = "polar";
                    else group = "nonPolar";
                    writer.print("\t<" + label + " name=\"" + type + "\" group=\"" + group + "\" backbone=\"" + backbone + "\" nNowhere=\""+nNowhereAtoms[ordinal]+"\">\n");
                    if (lastCNC != 0)
                        for (int iBinCNC = 0; iBinCNC <= lastCNC; iBinCNC++) {
                            writer.print("\t\t" + "<CNCbin>" + iBinCNC);
                            for (int jBinHBC = 0; jBinHBC <= lastHBC; jBinHBC++)
                                writer.print("<HBCbin>  " + jBinHBC/2.0 + " " + histograms2D[ordinal][iBinCNC][jBinHBC] + "</HBCbin>\t");
                            writer.print("</CNCbin>\n");
                        }
                    writer.print("\t</" + label+ ">\n");
                }
            }
            writer.println("</" + label + "_statistics>");
            writer.close();
        }
    }
}

//        public void analyzeProtein(Protein protein) {
//            EnergyCreator[] creator = {new SolvationCreator()};
//            try {
//                TotalEnergy energy = new TotalEnergy(protein, creator, commands, "");
//                energy.update();
//                energy.evaluateAtoms();
//                double maxEnergy = -10000;
//                double minEnergy = 10000;
//                Atom minEnergyAtom = null, maxEnergyAtom = null;
//                SolvationEnergy se = (SolvationEnergy) energy.getEnergyTerm(new SolvationEnergy());
//                double[] hbc = se.hbc();
//                double[] cnc = se.cnc();
//                AtomData[] atomData = new AtomData[protein.atoms().size()];
//                for (Atom atom : protein.atoms()) {
//                    if (atom.energy() < minEnergy) {
//                    minEnergy = atom.energy();
//                    minEnergyAtom = atom;
//                }
//                if (atom.energy() > maxEnergy) {
//                    maxEnergy = atom.energy();
//                    maxEnergyAtom = atom;
//                }
//            }
//            for (Atom atom : protein.atoms()) {
//                double e = atom.energy();
//                atom.setTemperatureFactor(99.98 * (e - minEnergy) / (maxEnergy-minEnergy));
//                atomData[atom.number()] = new AtomData(atom,hbc[atom.number()],cnc[atom.number()],e);
//            }
//            System.out.println("minEnegyAtom = " + minEnergy + " " + hbc[minEnergyAtom.number()] + " " + cnc[minEnergyAtom.number()] + " " + minEnergyAtom);
//            System.out.println("maxEnegyAtom = " + maxEnergy + " " + hbc[maxEnergyAtom.number()] + " " + cnc[maxEnergyAtom.number()] + " " + maxEnergyAtom);
//            MeshiWriter writer = new MeshiWriter(protein.name() + ".out.pdb");
//            protein.atoms().print(writer);
//            writer.close();
//
//            Arrays.sort(atomData);
//            writer = new MeshiWriter(protein.name()+".atomData.txt");
//            for (AtomData ad : atomData)
//                writer.println(ad);
//            writer.close();
//        } catch (Exception ex) {
//            ex.printStackTrace();
//            System.out.println("------------------------");
//            throw new RuntimeException(ex);
//        }
//    }

//    private static class AtomData implements Comparable<AtomData>{
//        Atom atom;
//        double hbc,cnc, energy;
//
//        AtomData(Atom atom, double hbc, double cnc, double energy){
//            this.atom = atom;
//            this.hbc = hbc;
//            this.cnc = cnc;
//            this.energy = energy;
//        }
//        public String toString() {
//            return "AtomData: "+atom+" "+cnc+" "+hbc+" "+energy;
//        }
//
//        public int compareTo(AtomData other) {
//            if (energy < other.energy) return -1;
//            if (energy > other.energy) return 1;
//            return 0;
//        }
//    }



//        public static void analyzeTarget(String outputDirectoryName) throws IOException, AlignmentException, UpdateableException, EvaluationException {
//        outputDir = new File(outputDirectoryName);
//        if (!outputDir.exists())
//            outputDir.mkdir();
////        File benchmark = new File("CASP8_9_benchmark");
//        File benchmark = new File("CASP10");
//        File[] targetDirs = benchmark.listFiles();
//        for (File targetDir : targetDirs) {
//            if (
//                    targetDir.isDirectory() && targetDir.getName().startsWith("T0")) {
//                String targetName = targetDir.getName();
//                File outputDir = new File(outputDirectoryName + "\\" + targetName);
//                if (outputDir.exists()) continue;
//                outputDir.mkdir();
//                String nativeFileName = targetDir.getPath() + "\\" + targetName + ".N.pdb";
//                Protein nativeStructure = new Protein(new AtomList(nativeFileName), new ResidueExtendedAtomsCreator());
//                File[] fileNames = targetDir.listFiles();
//                for (File file : fileNames) {
//                    if (file.getName().endsWith("pdb")) {
//                        System.out.println("\n\nNow analyzing " + file.getPath());
//                        ModelToAnalyze modelToAnalyze = new ModelToAnalyze(file.getPath(), solvationCreator, contactsCreator, "dssp", commands);
//                        ModelAnalyzer modelAnalyzer = new ModelAnalyzer(modelToAnalyze.model, nativeStructure, modelToAnalyze.model, modelToAnalyze.energy, ResidueAlignmentMethod.IDENTITY);
//                        MeshiInfoXMLwriter writer = new MeshiInfoXMLwriter(outputDirectoryName + "/" + targetName + "/" + file.getName() + ".xml");
//                        ProteinInfoOLd info = modelAnalyzer.analyze("Atom environment features");
//                        ProteinInfoListOld infoList = new ProteinInfoListOld(modelToAnalyze.model.name());
//                        infoList.add(info);
//                        infoList.print(writer);
//                        writer.close();
//                    }
//                }
//            }
//        }
//    }
//        double[] energies = modelToAnalyze.solvationEnergy.atomEnergies;
//        double[] atomTypeEnergies = new double[AtomType.values().length];
//        for (int i = 1; i< AtomType.values().length; i++)
//            atomTypeEnergies[i] = 0;
//
//        atomTypeEnergies[ordinal] += energies[atomCore.number];
//        modelsDir = new File(outputDir.getPath() + "/models");
//        modelNames = new MeshiWriter(modelsDir.getAbsolutePath() + "\\modelsNames.txt");
//        targetSolvationData = new MeshiWriter(modelsDir.getAbsolutePath() + "\\modelsEnergies.dat");
//        first = false;
//        targetSolvationData.print(modelNumber + "\t" + model.chain().getNonDummyResidues().size() + "\t" + gdt[0] + "   " + solvationInfoElement.value() + "\t" + "  " + solvationInfoElement.energyInfo.value() + "\t");
//        targetSolvationData.print(solvationInfoElement.bbPolarInfo.value() + "\t" + solvationInfoElement.bbCarbonInfo.value() + "\t" + solvationInfoElement.scPolarInfo.value() + "\t" + solvationInfoElement.scCarbonInfo.value()+"\t"+solvationInfoElement.hbInfo.value()+"\t"+solvationInfoElement.stdInfo.value()+"\t"+solvationInfoElement.entropyInfo.value()+"\t");
//        targetSolvationData.print(solvationInfoElement.bbPolarN.value() + "\t" + solvationInfoElement.bbCarbonN.value() + "\t" + solvationInfoElement.scPolarN.value() + "\t" + solvationInfoElement.scCarbonN.value());
//        targetSolvationData.print("\t"+contactsInfo.contacts12.value()+"\t"+contactsInfo.contacts14.value()+"\t"+contactsInfo.contacts14CoreRatio.value()+"\t"+contactsInfo.sasaRatio.value());





//    public static void removeOutliers(Protein protein) {
//        for (Atom atom:protein.atoms())  {
//            if ((!atom.nowhere()) && (atom.x() < -999))
//                atom.setNowhere();
//        }
//    }

//    public static void testProtein(Protein protein,SolvationEnergy solvationEnergy) throws IOException{
//        MeshiWriter writer = new MeshiWriter("proteinTest.dat");
//        solvationEnergy.evaluate();
//        double[] cnc = solvationEnergy.cnc();
//        double[] hbc = solvationEnergy.hbc();
//        double[] energies = solvationEnergy.atomEnergies;
//        MolecularSystem molecularSystem = protein.atoms().molecularSystem;
//        testAtom(writer, molecularSystem, cnc, hbc, solvationEnergy, 1);
//        for (int iAtom = 0; iAtom < molecularSystem.size(); iAtom++) {
//            writer.println(molecularSystem.get(iAtom).type()+", "+molecularSystem.get(iAtom).atom + "  " + cnc[iAtom] + "   " + hbc[iAtom] + "   " + energies[iAtom]);
//        }
//        writer.close();
//    }

//    static void testAtom(MeshiWriter writer, MolecularSystem molecularSystem, double[] CNC, double[] HBC, SolvationEnergy solvationEnergy, int atomNumber){
//        double cnc = CNC[atomNumber];
//        double hbc = HBC[atomNumber];
//        AtomCore atom = molecularSystem.get(atomNumber);
//        Spline2D spline2D = solvationEnergy.splines()[atomNumber];
//        writer.println("Testing "+cnc+" , "+hbc+" , "+spline2D+" , "+atom);
//        spline2D.calc(cnc,hbc);
//        writer.println("spline2D.s = "+spline2D.s);
//        spline2D.test(writer, cnc, hbc);
//    }

//    private static class ModelToAnalyze {
//        public Protein model;
//        TotalEnergy energy;
//        public ModelToAnalyze(String fileName,
//                              SolvationCreator solvationCreator,
//                              ContactsCreator contactsCreator,
//                              String dsspDir,
//                              CommandList commands) throws AlignmentException, UpdateableException, IOException, EvaluationException {
//            new MolecularSystem();
//            model = new Protein(new AtomList(fileName), new ResidueExtendedAtomsCreator());
//            for (Atom atom : model.atoms()) {
//                if((!atom.nowhere()) && (atom.x() < -999))
//                    throw new RuntimeException("This is weird.\n"+atom+"  "+atom.x()+" "+(atom.x() < -999));
//            }
//            Utils.addHydrogens(model,commands);
//            removeOutliers(model);
//            fileName = fileName.substring(fileName.indexOf("T0"));
//            if (!Utils.AssignDSSP(model, dsspDir+"/"+fileName+".dssp")) return;
//            for (int iResidue = 0; iResidue < model.residues().size(); iResidue++) {
//                Residue residue = model.residues().get(iResidue);
//                if (residue.getSecondaryStructure() == SecondaryStructure.UNK) {
//                    residue.setSecondaryStructure(SecondaryStructure.COIL);
//                }
//            }
//            EnergyCreator[] energyCreators = {solvationCreator,contactsCreator, goapCreator};
//            energy = new TotalEnergy(model,energyCreators, DistanceMatrix.DistanceMatrixType.STANDARD, commands,"High and low resolution solvation terms");
//            energy.update();
//            energy.evaluate();
//        }
//    }


