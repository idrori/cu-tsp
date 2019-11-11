/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

import meshi.energy.EvaluationException;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateCreator;
import meshi.energy.simpleEnergyTerms.inflate.inflate.InflateType;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.geometry.ArcCos;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.OptimizerException;
import meshi.parameters.MeshiPotential;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.Filter;

import java.io.File;
import java.io.IOException;



public class LocalStructureAnalyzer extends MeshiProgram implements KeyWords {


    public static final String NAME = "LocalStructureAnalyzer";
    private static String inFileName, nativeFileName, outFileName,dsspFile;
    private static Protein model;
    private static CommandList commands;
    private static int seed;
    private static InflateCreator inflateByOtherModelCreator = new InflateCreator(InflateType.BY_OTHER_MODEL, new File("."), new MyFileFilter());
    private static String parentString;


    public static void main(String[] argv) throws IOException, OptimizerException, UpdateableException, EvaluationException {
        init(argv);
        ArcCos.useFastArcCos();
        model = null;

        parentString = getParentString(inFileName);
        model = Utils.getProtein(commands, inFileName, ResidueExtendedAtomsCreator.creator, Utils.defaultExceptionHandler);
        Utils.alignToX(model);
        addAtoms(model, commands);
        model.resetBonds();


        Utils.println("nAtoms " + model.atoms().size());
        Utils.println("nResidues " + model.residues().size());
        ((MyFileFilter) inflateByOtherModelCreator.filter).reset(model);

//      Secondary structure;
        if (!dsspFile.equals("NONE")) {
            Utils.AssignFullDSSP(model, dsspFile);
        }


        /*DSSPFull dssp = new DSSPFull(dsspFile);
        System.out.println("DDDDDDDDD = "+dssp.getSSInOneLine());
        */
        System.out.println("Secondary structure:");
        Utils.setSS(model, commands);

        String seqResidues = "";
        String seqSS = "";
        String seqSS3 = "";
        String seqSS7 = "";
        String seqSTR2 = "";

        Chain residues = model.chain();
        int iResStart=residues.firstNonDummyResidue().number();
        for (int iRes=iResStart; iRes < residues.size(); iRes ++){
            Residue residue = residues.residueAt(iRes);
            seqResidues+=residue.type().nameOneLetter();
            seqSS+=residue.getSecondaryStructure().getNameOneLetter();
            seqSS3+=residue.getLocalStructure().getDSSP3Reduction().getNameOneLetter();
            seqSS7+=residue.getLocalStructure().getDSSP7Reduction().getNameOneLetter();
            seqSTR2+=residue.getLocalStructure().getSTR2Reduction().getNameOneLetter();
        }

        System.out.println(seqResidues);
        System.out.println(seqSS);
        System.out.println(seqSS3);
        System.out.println(seqSS7);
        System.out.println(seqSTR2);

        MeshiWriter result = new MeshiWriter(outFileName);
        result.println(seqResidues);
        result.println(seqSS);
        result.println(seqSS3);
        result.println(seqSS7);
        result.println(seqSTR2);
        result.close();
    }

    public static void addAtoms(Protein model, CommandList commands) throws IOException{
        Command command = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY);
        String parametersDirectory = command.secondWord();
        BondParametersList bondParametersList  = new BondParametersList(parametersDirectory+"/" + MeshiPotential.BOND_PARAMETERS);
        AngleParametersList angleParametersList = new AngleParametersList(parametersDirectory+"/" + MeshiPotential.ANGLE_PARAMETERS);
        PlaneParametersList planeParametersList = new PlaneParametersList(parametersDirectory+"/" + MeshiPotential.PLANE_PARAMETERS);

        boolean hydrogenFailure = false;
        for (Atom atom : model.atoms()) {
            if (atom.isHydrogen() && atom.nowhere()) {
                PutHposLog puthLog = PutHpos.pos(atom, bondParametersList, angleParametersList);
                if (puthLog == null ) hydrogenFailure = true;
            }
        }
        AtomFinding.finalFindAtoms(model,bondParametersList,angleParametersList,planeParametersList);

    }

    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"commands", "inFileName", "dsspFile", "nativeStructure", "outFile", "seed"};
        String[] arguments = getArguments(keys, argv);

        seed = Integer.parseInt(arguments[5]);
        System.out.println("seed " + seed);
        initRandom(seed);

        commands = new CommandList(arguments[0]);
        inFileName = arguments[1];
        dsspFile   = arguments[2];
        nativeFileName = arguments[3];
        outFileName = arguments[4];
        if (commands.keyExists("verbose")) Utils.verboseOn();
        else Utils.verboseOff();
    }

    private static String getParentString(String fileName) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String out;
        while ((out = reader.readLine()) != null) {
            if (out.startsWith("PARENT"))
                return out;
        }
        return "PARENT N/A";
    }

    private static class MyFileFilter implements Filter {
        private String prefix = null;
        Protein thisModel;

        public void reset(Protein thisModel) {
            this.thisModel = thisModel;
            int index = thisModel.name().indexOf('.');
            if (index != -1) prefix = thisModel.name().substring(0, index);
            else prefix = thisModel.name();

        }

        public boolean accept(Object obj) {
            if (prefix == null) throw new RuntimeException("prefix is " + prefix);
            File file = (File) obj;
            String path = file.getAbsolutePath();
            if (!path.endsWith("pdb")) return false;
            if (path.indexOf("out") == -1) return false;
            if (file.getName().startsWith(prefix)) {
                try {
                    double rms = Rms.rms(thisModel, Protein.getCAproteinFromApdbFile(file), ResidueAlignmentMethod.IDENTITY);
                    if (rms < 1) return false;
                } catch (Exception ex ) {
                    writeFailureXml(model, ex);
                    throw new RuntimeException(ex.getMessage());}
                return true;
            }
            return false;
        }
    }

    private static void writeFailureXml(Protein model, Exception exception)  {
        MeshiWriter writer;
        try {
            writer = new MeshiWriter(outFileName + ".xml");
        } catch (IOException ex) {throw  new RuntimeException("Cannot write failure XML file after exception:\n"+exception+"\n"+"Due to "+ex);}
        writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\" ?> ");
        writer.println("<ProteinInfoList name=\"Failure report for Protein: " + model.sourceFile() + "\">");
        writer.print("<ProteinInfo  value=\"MCM_END\" step=\"0\" ");
//        for (Score scoreFunction : scoreFunctions) {
//            writer.print(scoreFunction.toString()+"_weightedMedianScore=\"0.05\" "+scoreFunction.toString()+"_interdecile=\"0\" ");
//        }
        writer.println("time=\"0\" fileName=\""+inFileName+"\" >");
        writer.println("<Exception>\n"+ exception + "\n</Exception>" );
        writer.println("</ProteinInfo>\n" + "</ProteinInfoList>\n");
        writer.close();
    }
}

