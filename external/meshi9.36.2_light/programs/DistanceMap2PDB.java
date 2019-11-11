package programs;

import meshi.energy.EvaluationException;
import meshi.energy.cooperativeDistanceMapConstraints.CooperativeDistanceMapConstraintsCreator;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.distanceMapConstraints.DistanceMapConstraints;
import meshi.energy.simpleEnergyTerms.distanceMapConstraints.DistanceMapConstraintsCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.repulsion.RepulsionCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.torsionConstraints.TorsionConstraints;
import meshi.energy.simpleEnergyTerms.torsionConstraints.TorsionConstraintsCreator;
import meshi.geometry.*;
import meshi.geometry.putH.PutHpos;
import meshi.geometry.putH.PutHposLog;
import meshi.molecularElements.Chain;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomBuilder;
import meshi.molecularElements.loops.AtomFinding;
import meshi.optimizers.LBFGS;
import meshi.energy.EnergyCreator;
import meshi.energy.TotalEnergy;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.Optimizer;
import meshi.optimizers.SteepestDecent;
import meshi.parameters.MeshiPotential;
import meshi.parameters.ResidueType;
import meshi.sequences.MeshiSequence;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.info.InfoType;
import org.json.simple.JSONArray;
import org.json.simple.JSONAware;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.File;
import java.io.IOException;
import java.io.Reader;

import static meshi.molecularElements.atoms.AtomBuilder.*;

public class DistanceMap2PDB extends MeshiProgram implements KeyWords {
    private static TorsionList torsionList;
    private static int ID;
    private static boolean refineFlag;
    private static String usageMessage = "\nUsage: DistanceMap2PDB <data json file> <repulsion json file> <commands file> <refine | noRefinement> \n";
    public static void main(String[] args) throws IOException, ParseException, UpdateableException, EvaluationException{
        if (args.length < 4)
            throw new RuntimeException(usageMessage);

        if ("refine".startsWith(args[3]))
            refineFlag = true;
        else {
            if ("noRefinement".startsWith(args[3]))
                refineFlag = false;
            else throw new RuntimeException(usageMessage);
        }


        JSONObject data = getData(args[0]);
        JSONArray repulsions = getRepultions(args[1]);
        String proteinName = (String) ((JSONObject) data.get("metaData")).get("proteinName");
        double distorsionStd = (Double) ((JSONObject) data.get("metaData")).get("distortionStd");
        double distanceRange = (Double) ((JSONObject) data.get("metaData")).get("distanceRange");
        ID = ((Long) ((JSONObject) data.get("metaData")).get("ID")).intValue();
        if (args.length == 5)
            initRandom(new Integer(args[4]));
        else
            initRandom(ID);
        CommandList commands = new CommandList(args[2]);

        String sequence = (String) data.get("proteinSequence");
        String mode;
        JSONObject distances = (JSONObject) data.get("caDistances");
        mode = "ca";
        if (distances == null) {
            distances = (JSONObject) data.get("cbDistances");
            mode = "cb";
        }
        JSONObject torsions = (JSONObject) data.get("torsions");
        MeshiWriter start = new MeshiWriter("start."+ID);
        start.close();
        Protein protein = reconstructProtein(proteinName, sequence, distances, torsions, repulsions, mode, commands);
        if (refineFlag)
            writeStructure(protein,"."+ID+".reconstructed");
    }

    private static void writeStructure(Protein protein, String sufix) throws IOException{
        File outFile = new File(protein.name()+sufix+".pdb");
        MeshiWriter writer = new MeshiWriter(outFile);
        protein.chain().print(writer);
        writer.close();
    }

    private static JSONObject getData(String fileName) throws ParseException, IOException, UpdateableException{
        Reader reader = new MeshiLineReader(fileName);
        JSONParser parser = new JSONParser();
        return (JSONObject) parser.parse(reader);
    }

    private static JSONArray getRepultions(String fileName) throws ParseException, IOException, UpdateableException{
        Reader reader = new MeshiLineReader(fileName);
        JSONParser parser = new JSONParser();
        JSONObject jsonObject =  (JSONObject) parser.parse(reader);
        return (JSONArray) jsonObject.get("repultions");
    }

    private static Protein reconstructProtein(String name, String sequence,
                                              JSONObject distances, JSONObject torsions, JSONArray repulsions, String mode,
                                              CommandList commands) throws UpdateableException, IOException, EvaluationException {
        MeshiSequence meshiSequence = new MeshiSequence(sequence, "");
        Protein protein = new Protein(meshiSequence, name, ResidueExtendedAtoms.creator);
        protein.residues().print();
        reconstructFromExtended(protein, torsions, distances, repulsions, mode, commands);
        writeStructure(protein,".initial.reconstructed."+ID);
        if (refineFlag) {
            addAtoms(protein, torsionList, commands);
            EnergyCreator[] energyCreators1 = {
                    new BondCreator(),
                    new AngleCreator(),
                    new PlaneCreator(),
                    new OutOfPlaneCreator(),
                    new RamachandranSidechainEnergyCreator(),
                    new AtomicPairwisePMFSummaCreator(),
                    new TorsionConstraintsCreator(torsions),
                    new CooperativeDistanceMapConstraintsCreator(distances, mode)
            };
            relax(protein, energyCreators1, commands);
            EnergyCreator[] energyCreators2 = {
                    new BondCreator(),
                    new AngleCreator(),
                    new PlaneCreator(),
                    new OutOfPlaneCreator(),
                    new RamachandranSidechainEnergyCreator(),
                    new AtomicPairwisePMFSummaCreator(),
                    new TorsionConstraintsCreator(torsions),
                    new CooperativeDistanceMapConstraintsCreator(distances, mode)
            };
            relax(protein, energyCreators2, commands);
            EnergyCreator[] energyCreators3 = {
                    new BondCreator(),
                    new AngleCreator(),
                    new PlaneCreator(),
                    new OutOfPlaneCreator(),
                    new AtomicPairwisePMFSummaCreator(),
                    new RamachandranSidechainEnergyCreator(),
                    new TetherCreator(InfoType.TETHER)
            };
            relax(protein, energyCreators3, commands);
        }
        return protein;
    }

    private static Protein reconstructFromExtended(Protein protein,
                                                   JSONObject torsions, JSONObject distances, JSONArray repulsions, String mode,
                                                   CommandList commands) throws IOException, EvaluationException, UpdateableException{
        for (int i = 50; i < protein.chain().size(); i = i + 50) {
            backboneFromTorsions(protein, torsions, i);
            writeStructure(protein,".start."+i);
            fold(protein, distances, torsions, repulsions, commands, i, mode);
            writeStructure(protein,"end."+i);
        }
        backboneFromTorsions(protein, torsions, protein.chain().size()-1);
//        writeStructure(protein,"start.last");
        fold(protein, distances, torsions, repulsions, commands,protein.chain().size()-1, mode);
        return protein;
    }

    private static void backboneFromTorsions(Protein protein, JSONObject torsions, int last) throws IOException, EvaluationException, UpdateableException {
        Chain chain = protein.chain();
        Residue residue1 = chain.firstNonDummyResidue();
        if (residue1.ca().nowhere()) {
            residue1.amideN().setXYZ(0, 0, 0);
            residue1.ca().setXYZ(0, C_N_len, 0);
            AtomBuilder.buildByAngle(N_CA_C_ang, CA_C_len, residue1.ca(), residue1.amideN(), residue1.carbonylC());
            if (residue1.type != ResidueType.GLY)
                AtomBuilder.buildByTorsion(CB_OUT_OF_PLAIN, N_CA_CB_ang, CA_CB_len,
                        residue1.ca(), residue1.amideN(), residue1.carbonylC(), residue1.cb());
            residue1.cb().setXYZ(residue1.ca().x(), residue1.ca().y(), residue1.ca().z() + 1);
            residue1.carbonylO().setXYZ(residue1.carbonylC().x(), residue1.carbonylC().y(), residue1.carbonylC().z() + 1);
        }
        Residue prevResidue = null;
        double[] prevOmega = {Math.PI, 0}, prevPhi = {-150*Math.PI/180, 0}, prevPsi = {170*Math.PI/180,0};
        double[] phi = new double[2], psi = new double[2], omega = new double[2];
        for (Residue residue : chain) {
            if (!residue.dummy() & (residue.number() <= last)) {
                if (residue.ca().nowhere()) {
                    JSONArray residueTorsionsArray = (JSONArray) torsions.get("" + residue.number());
                    if (residueTorsionsArray != null) {
                        phi[0] = (double) ((JSONArray) residueTorsionsArray.get(0)).get(0);
                        psi[0] = (double) ((JSONArray) residueTorsionsArray.get(1)).get(0);
                        omega[0] = (double) ((JSONArray) residueTorsionsArray.get(2)).get(0);
                    } else {
                        phi[0] = -150 * Math.PI / 180;
                        phi[1] = 0;
                        psi[0] = 170 * Math.PI / 180;
                        psi[1] = 0;
                        omega[0] = Math.PI;
                        omega[1] = 0;
                    }
                    if ((prevResidue != null) && (!prevResidue.dummy()) & residue.amideN().nowhere()) {
                        //$N =gen_atom_LAT($topolines[$i-1][1]*PI/180,$CA_C_N_ang,$C_N_len,$pC,$pCA,$pN);
                        AtomBuilder.buildByTorsion(prevPsi[0], CA_C_N_ang, C_N_len,
                                prevResidue.carbonylC(), prevResidue.ca(), prevResidue.amideN(), residue.amideN());
                        //$CA=gen_atom_LAT(PI,$CA_C_N_ang,$N_CA_len,$N,$pC,$pCA);
                        AtomBuilder.buildByTorsion(prevOmega[0], CA_C_N_ang, N_CA_len,
                                residue.amideN(), prevResidue.carbonylC(), prevResidue.ca(), residue.ca());
                        //$C =gen_atom_LAT($topolines[$i][0]*PI/180,$N_CA_C_ang,$CA_C_len,$CA,$N,$pC);
                        AtomBuilder.buildByTorsion(phi[0], N_CA_C_ang, CA_C_len,
                                residue.ca(), residue.amideN(), prevResidue.carbonylC(), residue.carbonylC());
                        if (residue.type != ResidueType.GLY)
                            AtomBuilder.buildByTorsion(CB_OUT_OF_PLAIN, N_CA_CB_ang, CA_CB_len,
                                    residue.ca(), residue.amideN(), residue.carbonylC(), residue.cb());
                        AtomBuilder.buildByTorsion(-psi[0], CA_C_O_ang, C_O_len,
                                residue.ca(), residue.amideN(), residue.carbonylC(), residue.carbonylO());
                    }
                }
                prevResidue = residue;
                prevOmega = omega;
                prevPsi = psi;
            }
        }
    }

    private static void extendedBackbone(Protein protein, JSONObject torsions, JSONObject distances,CommandList commands) throws UpdateableException, IOException, EvaluationException {
        double x = -900;
        for (Residue residue : protein.chain()) {
            if (!residue.dummy()) {
                residue.ca().setXYZ(x, 0, 0);
                residue.amideN().setXYZ(x - 1.2, 0.6, randomNumberGenerator().nextDouble() / 10);
                residue.carbonylC().setXYZ(x + 1.2, 1.2, randomNumberGenerator().nextDouble() / 10);
                x += 3.5;
            }
        }
    }


    public static void addAtoms(Protein model, TorsionList torsionList, CommandList commands) throws IOException{
        Residue prevResidue = null;
        for (Residue residue : model.chain()) {
            if (!residue.dummy() & ((prevResidue != null) && (!prevResidue.dummy()))) {
                ResidueTorsions residueTorsions = new ResidueTorsions(residue, torsionList);
                Torsion phi = residueTorsions.phi();
                if (phi != null)
                    AtomBuilder.buildByTorsion(phi.torsion()-Math.PI, CA_C_O_ang, C_O_len,
                        residue.carbonylC(), residue.ca(), residue.amideN(), residue.carbonylO());
                else
                    AtomBuilder.buildByAngle(CA_C_O_ang, C_O_len, residue.carbonylC(), residue.ca(), residue.carbonylO());
                if (residue.type != ResidueType.PRO)
                    AtomBuilder.buildByTorsion(Math.PI, H_N_CA_ang, H_N_len,
                            residue.amideN(), residue.carbonylC(), residue.carbonylO(), residue.amideH());
                else
                    AtomBuilder.buildByTorsion(Math.PI, H_N_CA_ang, H_N_len,
                            residue.amideN(), residue.carbonylC(), residue.carbonylO(), residue.cd());
                if (residue.type != ResidueType.GLY )
                    AtomBuilder.buildByTorsion(CB_OUT_OF_PLAIN, N_CA_CB_ang, CA_CB_len,
                            residue.ca(), residue.amideN(), residue.carbonylC(), residue.cb());
            }
            prevResidue = residue;
        }
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


    private static void relax(Protein protein, EnergyCreator[] energyCreators, CommandList commands) throws UpdateableException{
        protein.atoms().get(0).molecularSystem.terminator().reset();
        TotalEnergy firstRelaxEnergy = new TotalEnergy(protein, energyCreators,
                DistanceMatrix.DistanceMatrixType.STANDARD, commands,"relaxEnergy");
        SteepestDecent steepestDecent = new SteepestDecent(firstRelaxEnergy,0.5,1000,200,0.00000001,0.5,1.5);
        try {
            steepestDecent.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        LBFGS lbfgs = Utils.getLBFGS(firstRelaxEnergy, commands, RELAX);
        try {
            lbfgs.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        Utils.println("nAtoms " + protein.atoms().size());
        Utils.println("nResidues " + protein.residues().size());//
    }

    public static void fold(Protein protein, JSONObject distances, JSONObject torsions, JSONArray repulsions,
                            CommandList commands, int last, String mode) throws UpdateableException, IOException {
        DistanceMapConstraintsCreator distanceMapConstraintsCreator = new DistanceMapConstraintsCreator(distances, mode);
        TorsionConstraintsCreator torsionConstraintsCreator = new TorsionConstraintsCreator(torsions);
        CooperativeDistanceMapConstraintsCreator cooperativeDistanceMapConstraintsCreator = new CooperativeDistanceMapConstraintsCreator(distances, mode);
        RepulsionCreator repulsionCreator = new RepulsionCreator(repulsions);
        EnergyCreator[] energyCreators = {
                new BondCreator("simple"),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                torsionConstraintsCreator,
                distanceMapConstraintsCreator,
                new RamachandranSidechainEnergyCreator(),
                repulsionCreator,
                cooperativeDistanceMapConstraintsCreator // Must be after distanceMapConstraintsCreator to ensure distances update.
        };
        int converges = 0;
        for (int i = 1; i <= 20 & (converges < 20); i++) {
            protein.atoms().get(0).molecularSystem.terminator().reset();
            TotalEnergy energy = new TotalEnergy(protein, energyCreators, commands, "foldEnergy");
            torsionList = ((TorsionConstraints) torsionConstraintsCreator.term()).torsionList();
            ((TorsionConstraints) torsionConstraintsCreator.term()).scaleWeights(2);
                Utils.println("---------------- Length: "+ last + "; Iteration " + i + " -----------------------------");
                ((DistanceMapConstraints) distanceMapConstraintsCreator.term()).setBeta(1 /  i);
                SteepestDecent steepestDecent = new SteepestDecent(energy, 0.5, 1000, 200, 0.00000001, 0.5, 1.5);
                try {
                    steepestDecent.run();
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
                LBFGS lbfgs = Utils.getLBFGS(energy, commands, MINIMIZE);
                try {
                    Optimizer.OptimizerStatus optimizerStatus = lbfgs.run();
                    if (optimizerStatus == Optimizer.OptimizerStatus.CONVERGED) {
                        converges++;
                        repulsionCreator.term().off();
                    }
                    else
                        converges = 0;
                } catch (Exception ex) {
                    throw new RuntimeException(ex);
                }
                Utils.println("nAtoms " + protein.atoms().size());
                Utils.println("nResidues " + protein.residues().size());
                Utils.alignToX(protein);
        }
    }

    //-------------------------------------- segments ---------------------------------------------------------
   /* private static void relaxSegments(Protein protein, JSONObject torsions, JSONObject distances, CommandList commands) throws UpdateableException{
        RepulsionCreator distanceMapConstraintsCreator = new RepulsionCreator(distances);
        EnergyCreator[] energyCreators = {
                new BondCreator("simple"),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                //     new RamachandranSidechainEnergyCreator(),
                new AtomicPairwisePMFSummaCreator(),
                new TorsionConstraintsCreator(torsions),
                distanceMapConstraintsCreator
        };
        TotalEnergy firstRelaxEnergy = new TotalEnergy(protein, energyCreators, commands,"relaxEnergy");
        ((DistanceMapConstraints) distanceMapConstraintsCreator.term()).setBeta(0.01);
        SteepestDecent steepestDecent = new SteepestDecent(firstRelaxEnergy,0.01,1000,200,0.00000001,0.5,1.5);
        try {
            steepestDecent.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        LBFGS lbfgs = Utils.getLBFGS(firstRelaxEnergy, commands, RELAX);
        try {
            lbfgs.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        Utils.println("nAtoms " + protein.atoms().size());
        Utils.println("nResidues " + protein.residues().size());//
    }

    private static Protein reconstructBySegments(Protein protein,
                                                 JSONObject torsions, JSONObject distances,
                                                 CommandList commands) throws IOException, EvaluationException, UpdateableException{
        Chain chain = protein.chain();
        int segmentLength = 10;
        for (int i = 1; i < chain.size(); i = i + segmentLength) {
            Utils.println("---------- Building segment " + i + " to " + (i + segmentLength - 1) + "  -------------------");
            buildSegment(protein, torsions, distances, commands, i, i + segmentLength - 1);
            MeshiWriter tmpWriter = new MeshiWriter("tmp.pdb");
            protein.chain().print(tmpWriter);
            tmpWriter.close();
        }
        fold(protein, distances, torsions, commands);
        return protein;
    }

    private static  void buildSegment(Protein protein,
                                      JSONObject torsions, JSONObject distances,
                                      CommandList commands, int first, int last) throws IOException, EvaluationException, UpdateableException {
        Random rnd = randomNumberGenerator();
        Chain chain = protein.chain();
        if (last >= chain.size())
            last = chain.size()-1;
        double x;
        if (first == 1)
            x = 0;
        else
            x = chain.residueAt(first - 1).ca().x() + 3.5;

        for (int iResidue = first; iResidue <= last; iResidue++) {
            Residue residue = chain.residueAt(iResidue);
            if (!residue.dummy()) {
                residue.ca().setXYZ(x, 0, 0);
                residue.amideN().setXYZ(x - 1.2, 0.6, rnd.nextDouble()/10);
                if (residue.amideH() != null)
                    residue.amideH().setXYZ(x - 1.2, -0.4, rnd.nextDouble()/10);
                residue.carbonylC().setXYZ(x + 1.2, 0.9, rnd.nextDouble()/10);
                residue.carbonylO().setXYZ(x + 1.2, 1.6, rnd.nextDouble()/10);
                x += 3.5;
            }
        }
        relaxSegments(protein, torsions, distances, commands);
    }
*/

}
