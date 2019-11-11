package programs;
import meshi.PDB.PdbReader;
import meshi.energy.simpleEnergyTerms.compositeTorsions.ResidueTorsions;
import meshi.geometry.*;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.ResidueType;
import meshi.sequences.*;
import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;
import org.json.simple.JSONObject;

import meshi.util.file.MeshiLineReader;

import javax.sound.midi.Sequence;
import java.io.IOException;
import java.util.Random;

public class Pdb2distanceMap extends MeshiProgram {
    private static boolean torsionFlag;
    private enum Mode {CA, CB};
    private static Mode mode;
    public static void main(String[] args) throws IOException{
        if (args.length != 5)
            throw new RuntimeException("Usage: Pdb2distanceMap <protein name> <distance range> <distortion std> <torsionsOn | torsionsOff>  <CA | CB>");
        initRandom();
        int ID = randomNumberGenerator().nextInt(10000);
        String proteinName = args[0];
        double distanceRange = new Double(args[1]);
        double distortionStd = new Double(args[2]);
        if (args[3].equals("torsionsOn"))
            torsionFlag = true;
        else
            torsionFlag = false;
        if (args[4].equals("CA"))
            mode = Mode.CA;
        else
            mode = Mode.CB;
        PdbReader pdbReader = new PdbReader(proteinName+".pdb");
        Protein protein = new Protein(new AtomList(pdbReader), ResidueExtendedAtoms.creator);
        MeshiSequence sequence = (new SequenceList(proteinName+".fasta.txt")).get(0);
        MeshiWriter writer = new MeshiWriter(proteinName+"._"+distanceRange+"_"+distortionStd+"_."+mode+"."+ID+".structData.json");
        JSONObject metaData = getMetaData(proteinName, distanceRange, distortionStd, ID);
        JSONObject jsonObject = new JSONObject();
        jsonObject.put("metaData",metaData);
        analyze(jsonObject, protein, sequence, distanceRange, distortionStd);
        writer.print(jsonObject);
        writer.close();
    }

    private static JSONObject getMetaData(String proteinName, double distanceRange, double distortionStd, int ID) {
        JSONObject jsonObject = new JSONObject();
        jsonObject.put("proteinName",proteinName);
        jsonObject.put("distanceRange", distanceRange);
        jsonObject.put("distortionStd", distortionStd);
        jsonObject.put("ID",ID);
        jsonObject.put("Mode",mode.toString());
        return jsonObject;
    }
    private static void analyze(JSONObject jsonObject, Protein protein, MeshiSequence sequence, double distanceRange, double distortionStd) {
        MolecularSystem molecularSystem = protein.molecularSystem;
        molecularSystem.terminator().reset();
        DistanceMatrix distanceMatrix = new DistanceMatrix(molecularSystem);
        jsonObject.put("proteinName",getName(protein));
        JSONObject sequenceObject = getSequence(protein, sequence);
        jsonObject.put("proteinSequence", sequenceObject.get("sequence"));
        if (mode == Mode.CA)
            jsonObject.put("caDistances",getDistances(protein, -1 * (int) sequenceObject.get("sequenceShift"), distanceRange, distortionStd));
        else
            jsonObject.put("cbDistances",getDistances(protein, -1 * (int) sequenceObject.get("sequenceShift"), distanceRange, distortionStd));
        jsonObject.put("torsions", getTorsions(protein, -1 * (int) sequenceObject.get("sequenceShift"), distanceMatrix, 0));
    }
    private static String getName(Protein protein) {
        return protein.name().substring(0,protein.name().length()-4);
    }

    private static JSONObject getSequence(Protein protein, MeshiSequence sequence) {
        JSONObject jsonObject = new JSONObject();
        SequenceAlignment alignment = SequenceAlignment.identityAlignment(protein.chain().sequence(), sequence);
        boolean found = false;
        int i = 0;
        SequenceAlignmentColumn column = null;
        while (!found) {
            column = alignment.get(i);
            AlignmentCell proteinCell = column.cell0();
            if ((!proteinCell.gap()) & (proteinCell.number >= 1))
                found = true;
        }
        int first = column.cell1().number() - column.cell0().number();

        String seqOut = "";
        for (AlignmentColumn clmn : sequence) {
            AlignmentCell sequenceCell = clmn.cell0();
            if (sequenceCell.number() - first >= 1)
                seqOut += sequenceCell.obj;
        }
        int sequenceShift;
        if (first < 0)
            sequenceShift = -first;
        else
            sequenceShift = 0;

        jsonObject.put("sequence",seqOut);
        jsonObject.put("sequenceShift",sequenceShift);
        return jsonObject;
    }
    private static JSONObject getDistances(Protein protein, int sequenceShift, double distanceRange, double distortionStd) {
        JSONObject distances = new JSONObject();
        Chain chain = protein.chain();
        for (int i = 1; i < chain.size(); i++) {
            Residue residueI = chain.residueAt(i);
            if (!residueI.dummy()) {
                Atom iAtom;
                if (mode == Mode.CA)
                    iAtom = residueI.ca();
                else {
                    if (residueI.cb() != null)
                        iAtom = residueI.cb();
                    else if (residueI.type == ResidueType.GLY)
                        iAtom = residueI.ca();
                    else throw new RuntimeException("This is weird");
                }
                if (! iAtom.nowhere()) {
                    double xi = iAtom.x();
                    double yi = iAtom.y();
                    double zi = iAtom.z();
                    JSONObject row = new JSONObject();
                    distances.put(i + sequenceShift, row);
                    for (int j = 1; j < chain.size(); j++) {
                        Residue residueJ = chain.residueAt(j);
                        if (!residueJ.dummy()) {
                            Atom jAtom;
                            if (mode == Mode.CA)
                                jAtom = residueJ.ca();
                            else {
                                if (residueJ.cb() != null)
                                    jAtom = residueJ.cb();
                                else if (residueJ.type == ResidueType.GLY)
                                    jAtom = residueJ.ca();
                                else throw new RuntimeException("This is weird");
                            }
                            if (! jAtom.nowhere()) {
                            double dx = jAtom.x() - xi;
                            double dy = jAtom.y() - yi;
                            double dz = jAtom.z() - zi;
                            double distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
                            if (Math.abs(residueI.number() - residueJ.number()) > 1) {
                                distance += distortion(distortionStd);
                            }
                            double[] distanceAndRange = {distance, distanceRange};
                            row.put(j + sequenceShift, distanceAndRange);
                        }
                    }
                }
                }
            }
        }
        return distances;
    }

    private static double distortion(double distortionStd) {
        return randomNumberGenerator().nextGaussian()*distortionStd;
    }


    private static JSONObject getTorsions(Protein protein, int sequenceShift, DistanceMatrix distanceMatrix, double torsionRange) {
        AtomPairList allBonds = new AtomPairList(protein.chains(),null);
        AtomPairList bonds = new AtomPairList();
        for (AtomPair bond : allBonds) {
            if ((!bond.atom1().nowhere()) & (!bond.atom2().nowhere()))
                bonds.add(bond);
        }
        QuickAndDirtyAngleList angleList = new QuickAndDirtyAngleList(bonds,distanceMatrix);
        QuickAndDirtyTorsionList torsions = new QuickAndDirtyTorsionList(angleList,distanceMatrix);
        JSONObject torsionsOut = new JSONObject();

        for (Residue residue : protein.chain()) {
            if (! residue.dummy() & torsionFlag) {
                ResidueTorsions residueTorsions = new ResidueTorsions(residue, torsions);
                Torsion phiTorsion = residueTorsions.phi();
                Torsion psiTorsion = residueTorsions.psi();
                Torsion omegaTorsion = residueTorsions.omega();
                double[] phi = new double[2], psi = new double[2], omega = new double[2];
                if (phiTorsion != null) {
                    phi[0] = phiTorsion.torsion();
                    phi[1] = torsionRange;
                }
                else {
                    phi[0] = -90*Math.PI/180;
                    phi[1] = 160*Math.PI/180;
                }
                if (psiTorsion != null) {
                    psi[0] = psiTorsion.torsion();
                    psi[1] = torsionRange;
                }
                else{
                        psi[0] = 50*Math.PI/180;
                        psi[1] = 240*Math.PI/180;
                    }
                if (omegaTorsion != null) {
                    omega[0] = omegaTorsion.torsion();
                    omega[1] = torsionRange;
                }
                else {
                    omega[0] = Math.PI / 180;
                    omega[1] = 0;
                }
                double[][] torsionsAndRanges = {phi, psi, omega};
                torsionsOut.put(residue.number()+sequenceShift, torsionsAndRanges);
            }
        }
        return torsionsOut;
    }
}


/*
class JsonEncodeDemo {

    public static void main(String[] args) {
        JSONObject obj = new JSONObject();

        obj.put("name", "foo");
        obj.put("num", new Integer(100));
        obj.put("balance", new Double(1000.21));
        obj.put("is_vip", new Boolean(true));

        System.out.print(obj);
    }
}
*/