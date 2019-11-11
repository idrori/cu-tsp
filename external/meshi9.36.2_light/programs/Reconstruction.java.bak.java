package programs;

import meshi.energy.EnergyCreator;
import meshi.energy.EvaluationException;
import meshi.energy.TotalEnergy;
import meshi.energy.pairwiseNonBondedTerms.atomicPairwisePMFSumma.AtomicPairwisePMFSummaCreator;
import meshi.energy.simpleEnergyTerms.angle.AngleCreator;
import meshi.energy.simpleEnergyTerms.bond.BondCreator;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ramachandranSidechain.RamachandranSidechainEnergyCreator;
import meshi.energy.simpleEnergyTerms.distanceMapConstraints.DistanceMapConstraintsCreator;
import meshi.energy.simpleEnergyTerms.outOfPlane.OutOfPlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneCreator;
import meshi.energy.simpleEnergyTerms.plane.PlaneParametersList;
import meshi.energy.simpleEnergyTerms.repulsion.RepulsionCreator;
import meshi.energy.simpleEnergyTerms.tether.TetherCreator;
import meshi.energy.simpleEnergyTerms.torsionConstraints.TorsionConstraintsCreator;
import meshi.energy.solvation.SolvationCreator;
import meshi.geometry.AlphaCoordinates;
import meshi.geometry.Angle;
import meshi.geometry.Coordinates;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomBuilder;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.optimizers.LBFGS;
import meshi.optimizers.OptimizerException;
import meshi.optimizers.SteepestDecent;
import meshi.parameters.ResidueType;
import meshi.sequences.AlignmentException;
import meshi.sequences.MeshiSequence;
import meshi.sequences.ResidueAlignment;
import meshi.sequences.ResidueAlignmentMethod;
import meshi.util.*;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.filters.CbAndGlyFilter;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Date;

import static meshi.molecularElements.atoms.AtomBuilder.*;
import static meshi.util.info.InfoType.TETHER;
import static programs.Reconstruction.KEYS.*;


public class Reconstruction extends MeshiProgram implements KeyWords {
    public enum MODE {CA, CB};
    public MODE mode = null;
    public CommandList commands;
    public Protein protein, nativeStructure;
    MeshiWriter writer;
    XYZobject xyzObject;
    BackboneGeomerty backboneGeomerty;
    double[][][] phiPsiArray;
    ResidueAlignment alignmentToNative;
    TetherCreator tetherCreator;


    public static void main(String[] args) throws IOException, UpdateableException, ParseException, AlignmentException, EvaluationException, OptimizerException {
        long startTime = (new Date()).getTime();

        Reconstruction reconstruction = new Reconstruction(args);
        if (reconstruction.commands.keyExists("verbose"))
            Utils.verboseOff();

        Protein protein = reconstruction.reconstructProtein();

        protein.atoms().print(reconstruction.writer);
        reconstruction.writer.close();
        long endTime = (new Date()).getTime();
        System.out.println("Runtime = " + ((endTime - startTime) / 60000.0) + " minutes.");
    }


    public Reconstruction(String[] args) throws IOException, ParseException {
        getArguments(args);
    }


    private Protein reconstructProtein() throws AlignmentException, ParseException, UpdateableException, IOException, EvaluationException, OptimizerException {
        tetherCreator = null;
        Protein protein = initProtein();
        report1(protein, "Initial Rms & GDT = ");

        fixBackbone(protein, backboneGeomerty, phiPsiArray, commands);
        addSidueChains(protein, commands);
        report2("Before relaxation Rms & GDT = ");
        relaxation(protein);
        report2("Final Rms & GDT = ");
        return protein;
    }

    private void relaxation(Protein protein)
            throws IOException, UpdateableException, EvaluationException, AlignmentException, OptimizerException {
        if (mode == MODE.CA)
            tetherCreator = new TetherCreator(TETHER);
        else if (mode == MODE.CB)
            tetherCreator.term().on();

        AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator = new AtomicPairwisePMFSummaCreator();
        RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator = new RamachandranSidechainEnergyCreator();
        EnergyCreator[] energyCreators = {
                new BondCreator(),
                new AngleCreator(),
                new PlaneCreator(),
                new OutOfPlaneCreator(),
                ramachandranSidechainEnergyCreator,
                atomicPairwisePMFSummaCreator,
                new SolvationCreator(),
                tetherCreator
        };
        TotalEnergy energy = new TotalEnergy(protein, energyCreators, commands, "backboneEnergy");
        SteepestDecent steepestDecent = new SteepestDecent(energy, 0.5, 1000, 200, 0.00000001, 0.5, 1.5);
        try {
            steepestDecent.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        LBFGS lbfgs = Utils.getLBFGS(energy, commands, MINIMIZE);
        lbfgs.run();
        tetherCreator.term().scaleWeight(40);
        ramachandranSidechainEnergyCreator.term().scaleWeight(0.5);
        atomicPairwisePMFSummaCreator.term().scaleWeight(2);
        lbfgs.run();
        tetherCreator.term().scaleWeight(0.001);
        lbfgs.run();
        if (nativeStructure != null) {
            double rms = Rms.rms(alignmentToNative, Rms.RmsType.CA);
            double[] gdt = Rms.gdt(alignmentToNative, nativeStructure.chain().numberOfNonDummyResidues());
            System.out.println("Final Rms & GDT = " + rms + " " + gdt[0]);
            writer.println(+rms + " " + gdt[0]);
        }
        writer.println(energy.reportHeader());
        writer.println(energy.report(lbfgs.maxSteps));
    }

    private void report1(Protein protein, String prefix) throws AlignmentException {
        if (nativeStructure != null)
            alignmentToNative = new ResidueAlignment(nativeStructure.chain(), "native",
                    protein.chain(), "model", ResidueAlignmentMethod.IDENTITY);
        else
            alignmentToNative = null;
        report2(prefix);
    }

    public void report2(String prefix) {
        if (alignmentToNative != null) {
            double rms = Rms.rms(alignmentToNative, Rms.RmsType.CA);
            double[] gdt = Rms.gdt(alignmentToNative, nativeStructure.chain().numberOfNonDummyResidues());
            System.out.println(prefix + rms + " " + gdt[0]);
            writer.println(prefix + rms + " " + gdt[0]);
        }
    }

    private Protein initProtein()
            throws UpdateableException, AlignmentException, OptimizerException, EvaluationException, IOException{
        MeshiSequence meshiSequence = new MeshiSequence(xyzObject.stringSequence, "");
        Protein protein = new Protein(meshiSequence, xyzObject.proteinName, ResidueExtendedAtoms.creator);
        Chain chain = protein.chain();
        if (mode == MODE.CA)
            initCA(protein, chain);
        else if (mode == MODE.CB){
            initCB(protein, chain);
            report1(protein, "After initCB Rms & GDT = ");
        }
        else
            throw new RuntimeException("Weird mode = "+mode+" may be either \"CA\" or \"CB\"");

        AlphaCoordinates[] alphaTorsions = getAlphaTorsions(chain);
        checkChirality(alphaTorsions, chain);
        phiPsiArray = initPhiPsiArray(alphaTorsions, backboneGeomerty, chain);
        buildInitialModel(alphaTorsions, phiPsiArray, chain);

        return protein;
    }

    private void initCA(Protein protein, Chain chain) {
        for (int i = 1; i < chain.size(); i++) {
            double[] coor = xyzObject.xyz.get(i - 1);
            chain.residueAt(i).ca().setXYZ(coor[0], coor[1], coor[2]);
        }
    }

    private void initCB(Protein protein, Chain chain)
            throws UpdateableException, AlignmentException, OptimizerException, EvaluationException, IOException{
        for (int i = 1; i < chain.size(); i++) {
            double[] coor = xyzObject.xyz.get(i - 1);
            if (chain.residueAt(i).type == ResidueType.GLY)
                chain.residueAt(i).ca().setXYZ(coor[0], coor[1], coor[2]);
            else {
                Coordinates cbCoordinates = new Coordinates(coor[0], coor[1], coor[2]);
                chain.residueAt(i).cb().setXYZ(cbCoordinates);
                chain.residueAt(i).ca().setXYZ(new Coordinates(cbCoordinates, AtomBuilder.CA_CB_len));
            }
        }
        report1(protein, "Initial Rms & GDT = ");
        tetherCreator = new TetherCreator(TETHER,CbAndGlyFilter.filter);
        EnergyCreator[] energyCreators = {
                new DistanceMapConstraintsCreator(chain, "cb"),
                tetherCreator
        };
        TotalEnergy energy = new TotalEnergy(protein, energyCreators, commands, "CbEnergy");
        SteepestDecent steepestDecent = new SteepestDecent(energy, 0.5, 1000, 200, 0.00000001, 0.5, 1.5);
        steepestDecent.run();
        tetherCreator.term().scaleWeight(100);
        LBFGS lbfgs = Utils.getLBFGS(energy, commands, MINIMIZE);
        lbfgs.run();
        tetherCreator.term().scaleWeight(100);
        lbfgs.run();
    }


    private static void checkChirality(AlphaCoordinates[] alphaTorsions, Chain chain) {
        double alphaEnergy = getAlphaEnergy(alphaTorsions, chain);
        if (alphaEnergy > 0.3) {
            mirror(chain);
            alphaTorsions = getAlphaTorsions(chain);
            Utils.println("Invering a D protein.");
        }
    }

    private static double[][][] initPhiPsiArray(AlphaCoordinates[] alphaTorsions, BackboneGeomerty backboneGeomerty, Chain chain) {
        double[][][] phiPsiArray = new double[alphaTorsions.length][2][2];
        for (int i = 0; i < phiPsiArray.length; i++) {
            phiPsiArray[i][0][0] = phiPsiArray[i][1][0] = -9999;
        }
        for (int i = 2; i < chain.size(); i++) {
            if (alphaTorsions[i].isValid() & alphaTorsions[i - 1].isValid()) {
                int[] indices = new int[5];
                indices[0] = (int) Math.floor((Angle.rad2deg(alphaTorsions[i - 1].getAlpha()) + 180) / 10);
                indices[1] = (int) Math.floor((Angle.rad2deg(alphaTorsions[i].getAlpha()) + 180) / 10);
                indices[2] = (int) Math.floor(Angle.rad2deg(alphaTorsions[i - 1].getAngle1()) / 10 - 8);
                indices[3] = (int) Math.floor(Angle.rad2deg(alphaTorsions[i - 1].getAngle0()) / 90);
                indices[4] = (int) Math.floor(Angle.rad2deg(alphaTorsions[i].getAngle1()) / 90);
                double phi = getValue(backboneGeomerty.phiPsi, indices, 0, 0);
                double psi = getValue(backboneGeomerty.phiPsi, indices, 1, 0);
                double phiStd = getValue(backboneGeomerty.phiPsi, indices, 0, 1);
                double psiStd = getValue(backboneGeomerty.phiPsi, indices, 1, 1);
                phiPsiArray[i][0][0] = phi;
                phiPsiArray[i][1][0] = psi;
                phiPsiArray[i][0][1] = phiStd;
                phiPsiArray[i][1][1] = psiStd;
            }
        }
        return phiPsiArray;
    }

    private static void buildInitialModel(AlphaCoordinates[] alphaTorsions, double[][][] phiPsiArray, Chain chain) {
        Residue residue, prevResidue, nextResidue;
        residue = chain.residueAt(1);
        residue.amideN().setXYZ(residue.ca().x() - 1, residue.ca().y(), residue.ca().z());
        residue = chain.residueAt(chain.size() - 1);
        residue.carbonylC().setXYZ(residue.ca().x() + 1, residue.ca().y(), residue.ca().z());


        for (int i = 2; i < chain.size(); i++) {
            residue = chain.residueAt(i);
            prevResidue = chain.residueAt(i - 1);
            double[] unitV = {residue.ca().x() - prevResidue.ca().x(), residue.ca().y() - prevResidue.ca().y(), residue.ca().z() - prevResidue.ca().z(),};
            double l = Math.sqrt(unitV[0] * unitV[0] + unitV[1] * unitV[1] + unitV[2] * unitV[2]);
            unitV[0] /= l;
            unitV[1] /= l;
            unitV[2] /= l;
            if (prevResidue.carbonylC().nowhere())
                prevResidue.carbonylC().setXYZ(prevResidue.ca().x() + CA_C_len * unitV[0],
                        prevResidue.ca().y() + CA_C_len * unitV[1],
                        prevResidue.ca().z() + CA_C_len * unitV[2]);
            if (residue.amideN().nowhere())
                residue.amideN().setXYZ(residue.ca().x() - N_CA_len * unitV[0] + 0.3 * randomNumberGenerator().nextDouble(),
                        residue.ca().y() - N_CA_len * unitV[1],
                        residue.ca().z() - N_CA_len * unitV[2]);

            if (alphaTorsions[i].isValid() & alphaTorsions[i - 1].isValid()) {
                nextResidue = chain.residueAt(i + 1);
                double phi = Math.PI * phiPsiArray[i][0][0] / 180;
                double psi = Math.PI * phiPsiArray[i][1][0] / 180;
                AtomBuilder.buildByTorsion(phi, CA_C_N_ang, CA_C_len, residue.ca(), residue.amideN(), prevResidue.carbonylC(), residue.carbonylC());
                AtomBuilder.buildByTorsion(psi, N_CA_C_ang, C_N_len, residue.carbonylC(), residue.ca(), residue.amideN(), nextResidue.amideN());

            }
        }
        for (int i = 1; i < chain.size(); i++) {
            residue = chain.residueAt(i);
            if ((residue.cb() != null) & (residue.cb() != null))  {
                AtomBuilder.buildByTorsion(CB_OUT_OF_PLAIN, N_CA_CB_ang, CA_CB_len,
                        residue.ca(), residue.amideN(), residue.carbonylC(), residue.cb());
            }
            residue.carbonylO().setXYZ(residue.carbonylC().x() + 1, residue.carbonylC().y(), residue.carbonylC().z());
            if (residue.amideH() != null)
                residue.amideH().setXYZ(residue.amideN().x() + 1, residue.carbonylC().y(), residue.carbonylC().z());
            else if (residue.cd() != null)
                residue.cd().setXYZ(residue.amideN().x() + 1, residue.carbonylC().y(), residue.carbonylC().z());
        }

    }

    public static void addSidueChains(Protein protein, CommandList commands) throws UpdateableException, EvaluationException, AlignmentException {
        PlaneCreator planeCreator = new PlaneCreator(commands);
        PlaneParametersList planeParametersList = (PlaneParametersList) planeCreator.parametersList();
        if (planeParametersList == null)
            throw new RuntimeException("This is weird.");
        Scmod.scmod(commands, protein, 1);
        Utils.addAtoms(protein, false, commands, planeParametersList, false);
        Scmod.scmod(commands, protein, 5);
    }

    public void fixBackbone(Protein protein, BackboneGeomerty backboneGeomerty, double[][][] phiPsiArray, CommandList commands)
            throws ParseException, UpdateableException, OptimizerException, EvaluationException, AlignmentException, IOException {
        if (mode == MODE.CA)
            tetherCreator = new TetherCreator(TETHER);
        PlaneCreator planeCreator = new PlaneCreator();
        BondCreator bondCreator = new BondCreator("simple");
        RepulsionCreator repulsionCreator = new RepulsionCreator(backboneGeomerty.repulsions);
        AtomicPairwisePMFSummaCreator atomicPairwisePMFSummaCreator = new AtomicPairwisePMFSummaCreator();
        TorsionConstraintsCreator torsionConstraintsCreator = new TorsionConstraintsCreator(phiPsiArray);
        RamachandranSidechainEnergyCreator ramachandranSidechainEnergyCreator = new RamachandranSidechainEnergyCreator();


        EnergyCreator[] energyCreators = {
                bondCreator,
                new AngleCreator(),
                planeCreator,
                new OutOfPlaneCreator(),
                repulsionCreator,
                torsionConstraintsCreator,
                ramachandranSidechainEnergyCreator,
                atomicPairwisePMFSummaCreator,
                tetherCreator,
        };
        TotalEnergy energy = new TotalEnergy(protein, energyCreators, commands, "backboneEnergy");
        if (mode == MODE.CB) {
            tetherCreator.term().on();
            tetherCreator.term().scaleWeight(0.0001);
        }
        SteepestDecent steepestDecent = new SteepestDecent(energy, 0.5, 1000, 200, 0.00000001, 0.5, 1.5);
        try {
            steepestDecent.run();
        } catch (Exception ex) {
            throw new RuntimeException(ex);
        }
        LBFGS lbfgs = Utils.getLBFGS(energy, commands, MINIMIZE);
        lbfgs.run();
        torsionConstraintsCreator.term().scaleWeight(2);
        tetherCreator.term().scaleWeight(40);
        bondCreator.term().scaleWeight(10);
        planeCreator.term().scaleWeight(10);
        //repulsionCreator.term().scaleWeight(0.1);
        ramachandranSidechainEnergyCreator.term().scaleWeight(1);
        atomicPairwisePMFSummaCreator.term().scaleWeight(1);
        lbfgs.run();
        ramachandranSidechainEnergyCreator.term().scaleWeight(1);
        atomicPairwisePMFSummaCreator.term().scaleWeight(1);
        tetherCreator.term().scaleWeight(20);
        bondCreator.term().scaleWeight(10);
        repulsionCreator.term().scaleWeight(0.0001);
        torsionConstraintsCreator.term().scaleWeight(0.5);
        lbfgs.run();
        torsionConstraintsCreator.term().scaleWeight(0.5);
        lbfgs.run();
        torsionConstraintsCreator.term().scaleWeight(0.5);
        lbfgs.run();
        torsionConstraintsCreator.term().scaleWeight(0.5);
        lbfgs.run();
        torsionConstraintsCreator.term().scaleWeight(0.5);
        lbfgs.run();
        torsionConstraintsCreator.term().scaleWeight(0.5);
        lbfgs.run();
        tetherCreator.term().scaleWeight(0.5);
        lbfgs.run();
        tetherCreator.term().scaleWeight(0.5);
        lbfgs.run();
        tetherCreator.term().scaleWeight(0.5);
        lbfgs.run();

    }

    private static double getValue(JSONArray phiPsi, int[] indices, int iPhiPsi, int iValueStd) {
        for (int i = 0; i < indices.length; i++) {
            if (indices[i] >= phiPsi.size())
                indices[i] = phiPsi.size() - 1;
            else if (indices[i] < 0)
                indices[i] = 0;
            phiPsi = (JSONArray) phiPsi.get(indices[i]);
        }
        phiPsi = (JSONArray) phiPsi.get(iPhiPsi);
        return (Double) phiPsi.get(iValueStd);
    }

    private static AlphaCoordinates[] getAlphaTorsions(Chain chain) {
        AlphaCoordinates[] alphaTorsions = new AlphaCoordinates[chain.size()];
        for (int i = 0; i < alphaTorsions.length; i++) {
            alphaTorsions[i] = new AlphaCoordinates();
        }
        chain.firstNonDummyResidue().ca().molecularSystem.terminator().reset();
        DistanceMatrix distanceMatrix = new DistanceMatrix(chain.firstNonDummyResidue().ca().molecularSystem);
        for (int i = 2; i < alphaTorsions.length - 2; i++) {
            int im1 = i - 1;
            int ip1 = i + 1;
            int ip2 = i + 2;
            Residue residueIm1 = chain.residueAt(im1);
            Residue residueI = chain.residueAt(i);
            Residue residueIp1 = chain.residueAt(ip1);
            Residue residueIp2 = chain.residueAt(ip2);
            alphaTorsions[i] = new AlphaCoordinates(residueIm1.ca(), residueI.ca(), residueIp1.ca(), residueIp2.ca(), distanceMatrix);
        }
        return alphaTorsions;
    }

    private static double getAlphaEnergy(AlphaCoordinates[] alphaTorsions, Chain chain) {
        double[] ENERGY = {-0.10, -0.08, -0.22, -0.26, -0.38, -0.45, -0.41, -0.50, -0.66, -0.74, -0.85, -0.89, -0.99,
                -0.98, -1.02, -1.03, -1.15, -1.13, -1.21, -1.27, -1.34, -1.22, -1.30, -1.36, -1.40, -1.34, -1.41, -1.50, -1.43, -1.42,
                -1.36, -1.65, -1.39, -1.55, -1.49, -1.42, -1.41, -1.51, -1.45, -1.42, -1.52, -1.32, -1.42, -1.44, -1.46, -1.60, -1.52,
                -1.35, -1.44, -1.46, -1.35, -1.43, -1.33, -1.23, -1.28, -1.31, -1.23, -1.30, -1.17, -1.22, -1.16, -1.04, -1.04, -1.22,
                -1.02, -0.99, -1.00, -1.03, -1.04, -0.96, -0.98, -1.08, -0.99, -0.92, -0.89, -0.87, -0.90, -0.81, -0.78, -0.78, -0.69,
                -0.65, -0.59, -0.53, -0.53, -0.35, -0.45, -0.45, -0.45, -0.33, -0.33, -0.20, -0.02, -0.06, 0.06, -0.02, 0.05, 0.32,
                0.28, 0.29, 0.27, 0.59, 0.54, 0.67, 0.82, 0.79, 1.05, 1.22, 1.30, 1.41, 1.44, 1.61, 1.69, 1.85, 1.95,
                2.17, 2.33, 2.22, 2.56, 2.62, 2.38, 2.68, 2.97, 3.21, 3.34, 3.60, 3.83, 4.19, 4.34, 4.44, 4.48, 4.44,
                4.44, 4.19, 4.33, 3.76, 3.55, 3.77, 3.20, 2.83, 3.02, 2.79, 2.73, 2.28, 2.44, 2.75, 2.40, 2.29, 1.99,
                2.22, 1.87, 1.97, 1.94, 1.73, 1.76, 1.64, 1.51, 1.55, 1.45, 1.39, 1.32, 1.15, 1.32, 1.12, 1.03, 0.86,
                0.83, 0.74, 0.83, 0.56, 0.41, 0.54, 0.44, 0.46, 0.31, 0.10, 0.27, 0.10, -0.01, -0.20, 0.20, 0.01, -0.10,
                -0.27, -0.10, -0.31, -0.46, -0.44, -0.54, -0.41, -0.56, -0.83, -0.74, -0.83, -0.86, -1.03, -1.12, -1.32, -1.15, -1.32,
                -1.39, -1.45, -1.55, -1.51, -1.64, -1.76, -1.73, -1.94, -1.97, -1.87, -2.22, -1.99, -2.29, -2.40, -2.75, -2.44, -2.28,
                -2.73, -2.79, -3.02, -2.83, -3.20, -3.77, -3.55, -3.76, -4.33, -4.19, -4.44, -4.44, -4.48, -4.44, -4.34, -4.19, -3.83,
                -3.60, -3.34, -3.21, -2.97, -2.68, -2.38, -2.62, -2.56, -2.22, -2.33, -2.17, -1.95, -1.85, -1.69, -1.61, -1.44, -1.41,
                -1.30, -1.22, -1.05, -0.79, -0.82, -0.67, -0.54, -0.59, -0.27, -0.29, -0.28, -0.32, -0.05, 0.02, -0.06, 0.06, 0.02,
                0.20, 0.33, 0.33, 0.45, 0.45, 0.45, 0.35, 0.53, 0.53, 0.59, 0.65, 0.69, 0.78, 0.78, 0.81, 0.90, 0.87,
                0.89, 0.92, 0.99, 1.08, 0.98, 0.96, 1.04, 1.03, 1.00, 0.99, 1.02, 1.22, 1.04, 1.04, 1.16, 1.22, 1.17,
                1.30, 1.23, 1.31, 1.28, 1.23, 1.33, 1.43, 1.35, 1.46, 1.44, 1.35, 1.52, 1.60, 1.46, 1.44, 1.42, 1.32,
                1.52, 1.42, 1.45, 1.51, 1.41, 1.42, 1.49, 1.55, 1.39, 1.65, 1.36, 1.42, 1.43, 1.50, 1.41, 1.34, 1.40,
                1.36, 1.30, 1.22, 1.34, 1.27, 1.21, 1.13, 1.15, 1.03, 1.02, 0.98, 0.99, 0.89, 0.85, 0.74, 0.66, 0.50,
                0.41, 0.45, 0.38, 0.26, 0.22, 0.08, 0.10};

        double energy = 0;
        int n = 0;
        for (int i = 1; i < chain.size(); i++) {
            Residue residue = chain.residueAt(i);
            if ((residue.type() != ResidueType.GLY) & (residue.type() != ResidueType.PRO) &
                    alphaTorsions[i].isValid()) {
                int index = (int) Math.floor(Angle.rad2deg(alphaTorsions[i].getAlpha()) + 180);
                energy += ENERGY[index];
                n++;
            }
        }
        return energy / n;
    }

    private static void mirror(Chain chain) {
        Utils.println("D coordinates - applying mirror symmetry.");
        for (Residue residue : chain) {
            for (Atom atom : residue.getAtoms()) {
                if (!atom.nowhere())
                    atom.setXYZ(-atom.x(), atom.y(), atom.z());
            }
        }
    }

    private static class XYZobject {
        String stringSequence;
        ArrayList<double[]> xyz = new ArrayList<>();
        String proteinName;

        public XYZobject(String xyzFileName) throws IOException {
            MeshiLineReader reader = new MeshiLineReader(xyzFileName);
            reader.readLine();
            proteinName = reader.readLine();
            String line = "";
            stringSequence = "";
            while (line != null) {
                line = reader.readLine();
                System.out.println(line);
                if (line != null) {
                    String[] words = line.split("\t");
                    if (words.length > 1) {
                        //System.out.println(words[0] + "," + words[1] + "," + words[2]+"," + words[3]);
                        stringSequence += words[0];
                        double[] coodinates = {(new Double(words[1])), (new Double(words[2])), (new Double(words[3]))};
                        xyz.add(coodinates);
                    }
                }
            }
        }
    }

    private static class BackboneGeomerty {
        JSONArray repulsions;
        JSONArray phiPsi;

        public BackboneGeomerty(String backBoneGeometryFileName) throws IOException, ParseException {
            Reader jsonReader = new MeshiLineReader(backBoneGeometryFileName);
            JSONParser parser = new JSONParser();
            JSONObject jsonObject = (JSONObject) parser.parse(jsonReader);
            repulsions = (JSONArray) jsonObject.get("minLocal");
            phiPsi = (JSONArray) jsonObject.get("phiPsi");

        }
    }

    private static void usageMessage(String message) {
        System.out.println(message);
        System.out.println("Usage: java Reconstruction <-key1=value1> <-key2=value2> ....");
        System.out.print("key list: ");
        for (KEYS key : KEYS.values())
            System.out.print(key + "  ");
        System.out.println();
        throw new RuntimeException("Cannot continue");
    }

//    private static class Arguments extends EnumMap<KEYS,Object>{
//        private enum KEYS {commands, xyzFile, nativeStructure, outFile, mode, seed, protocol};


    public enum KEYS {commands, xyzFile, nativeStructure, outFile, mode, seed, protocol, geometry, alphaTorsions, phiPsiArray, alignmentToNative}
    public void getArguments(String[] args) throws IOException, ParseException {
            String commandsFile = getArgument(KEYS.commands.toString(), args);
            if (commandsFile == null)
                usageMessage("Argument " + commands + " is mandatory yet missing.");
            commands = new CommandList(commandsFile);

            String xyz = getArgument(xyzFile.toString(), args);
            if (xyz == null)
                usageMessage("Argument " + xyzFile + " is mandatory yet missing.");
            xyzObject = new XYZobject(xyz);

            String nativeFileName = getArgument(KEYS.nativeStructure.toString(), args);
            if (nativeFileName != null)
                nativeStructure = new Protein(new AtomList(nativeFileName), ResidueExtendedAtoms.creator);
            else
                nativeStructure = null;

            String seedString = getArgument(KEYS.seed.toString(), args);
            if (seedString == null)
                usageMessage("Argument " + KEYS.seed + " is mandatory yet missing.");
            initRandom(new Integer(seedString));

            String outFileName = getArgument(outFile.toString(), args);
            if (outFileName == null)
                usageMessage("Argument " + outFile + " is mandatory yet missing.");
            writer = new MeshiWriter(outFileName + "." + seed() + ".pdb");

            String geometryFileName = getArgument(geometry.toString(), args);
            if (geometryFileName == null)
                usageMessage("Argument " + geometry + " is mandatory yet missing.");
            backboneGeomerty = new BackboneGeomerty(geometryFileName);

            String modeString = getArgument("mode", args);
            if (modeString == null)
                usageMessage("Argument "+mode+" is mandatory, yet missing.");
            if (modeString.equals("CA"))
                mode = MODE.CA;
            else if (modeString.equals("CB"))
                mode = MODE.CB;
            else
                throw new RuntimeException("Weird mode = "+mode+" may be either \"CA\" or \"CB\"");

        }
}
