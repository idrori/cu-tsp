package programs;

import meshi.PDB.PdbReader;
import meshi.energy.EnergyInfoElement;
import meshi.geometry.*;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.ResidueType;
import meshi.util.ProteinAnalyzer;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;
import meshi.util.info.MeshiInfo;
import meshi.util.info.ProteinInfo;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import static programs.ContactConstraints.ContactConstraintsAnalyzer.ENERGY;

/**
 */
public class ContactConstraints {

    public static void main(String[] args) throws Exception{
        ContactConstraintsAnalyzer analyzer = new ContactConstraintsAnalyzer();
        JSONObject residues = null;
        Utils.verboseOff();
        File[] files = (new File(".")).listFiles();
        int i = 1;
        for (File file : files) {
            PdbReader pdbReader = new PdbReader(file.getAbsolutePath());
            new MolecularSystem();
            System.out.println(file.getName());
            Protein protein = new Protein(new AtomList(pdbReader), ResidueExtendedAtoms.creator);
            if ( i == 2) {
                residues = getResidues(protein);
                analyzer.addResidues(residues);
            }
            analyzer.analyze(protein);
            analyzer.energyStatisticsE.doneProtein(protein);
            if (i % 200 == 0)
                analyzer.toJason(i+".json");
            i++;
        }
        analyzer.addResidues(residues);
        analyzer.toJason("contactConstraints.json");
    }

    static class ContactConstraintsAnalyzer implements ProteinAnalyzer {
        protected static double[] ENERGY =   {-0.10, -0.08, -0.22, -0.26, -0.38, -0.45, -0.41, -0.50, -0.66, -0.74, -0.85, -0.89, -0.99,
                -0.98, -1.02, -1.03, -1.15, -1.13, -1.21, -1.27, -1.34, -1.22, -1.30, -1.36, -1.40, -1.34, -1.41, -1.50, -1.43, -1.42,
                -1.36, -1.65, -1.39, -1.55, -1.49, -1.42, -1.41, -1.51, -1.45, -1.42, -1.52, -1.32, -1.42, -1.44, -1.46, -1.60, -1.52,
                -1.35, -1.44, -1.46, -1.35, -1.43, -1.33, -1.23, -1.28, -1.31, -1.23, -1.30, -1.17, -1.22, -1.16, -1.04, -1.04, -1.22,
                -1.02, -0.99, -1.00, -1.03, -1.04, -0.96, -0.98, -1.08, -0.99, -0.92, -0.89, -0.87, -0.90, -0.81, -0.78, -0.78, -0.69,
                -0.65, -0.59, -0.53, -0.53, -0.35, -0.45, -0.45, -0.45, -0.33, -0.33, -0.20, -0.02, -0.06,  0.06, -0.02,  0.05,  0.32,
                0.28,  0.29,  0.27,  0.59,  0.54,  0.67,  0.82,  0.79,  1.05,  1.22,  1.30,  1.41,  1.44,  1.61,  1.69,  1.85,  1.95,
                2.17,  2.33,  2.22,  2.56,  2.62,  2.38,  2.68,  2.97,  3.21,  3.34,  3.60,  3.83,  4.19,  4.34,  4.44,  4.48,  4.44,
                4.44,  4.19,  4.33,  3.76,  3.55,  3.77,  3.20,  2.83,  3.02,  2.79,  2.73,  2.28,  2.44,  2.75,  2.40,  2.29,  1.99,
                2.22,  1.87,  1.97,  1.94,  1.73,  1.76,  1.64,  1.51,  1.55,  1.45,  1.39,  1.32,  1.15,  1.32,  1.12,  1.03,  0.86,
                0.83,  0.74,  0.83,  0.56,  0.41,  0.54,  0.44,  0.46,  0.31,  0.10,  0.27,  0.10, -0.01, -0.20,  0.20,  0.01, -0.10,
                -0.27, -0.10, -0.31, -0.46, -0.44, -0.54, -0.41, -0.56, -0.83, -0.74, -0.83, -0.86, -1.03, -1.12, -1.32, -1.15, -1.32,
                -1.39, -1.45, -1.55, -1.51, -1.64, -1.76, -1.73, -1.94, -1.97, -1.87, -2.22, -1.99, -2.29, -2.40, -2.75, -2.44, -2.28,
                -2.73, -2.79, -3.02, -2.83, -3.20, -3.77, -3.55, -3.76, -4.33, -4.19, -4.44, -4.44, -4.48, -4.44, -4.34, -4.19, -3.83,
                -3.60, -3.34, -3.21, -2.97, -2.68, -2.38, -2.62, -2.56, -2.22, -2.33, -2.17, -1.95, -1.85, -1.69, -1.61, -1.44, -1.41,
                -1.30, -1.22, -1.05, -0.79, -0.82, -0.67, -0.54, -0.59, -0.27, -0.29, -0.28, -0.32, -0.05,  0.02, -0.06,  0.06,  0.02,
                0.20,  0.33,  0.33,  0.45,  0.45,  0.45,  0.35,  0.53,  0.53,  0.59,  0.65,  0.69,  0.78,  0.78,  0.81,  0.90,  0.87,
                0.89,  0.92,  0.99,  1.08,  0.98,  0.96,  1.04,  1.03,  1.00,  0.99,  1.02,  1.22,  1.04,  1.04,  1.16,  1.22,  1.17,
                1.30,  1.23,  1.31,  1.28,  1.23,  1.33,  1.43,  1.35,  1.46,  1.44,  1.35,  1.52,  1.60,  1.46,  1.44,  1.42,  1.32,
                1.52,  1.42,  1.45,  1.51,  1.41,  1.42,  1.49,  1.55,  1.39,  1.65,  1.36,  1.42,  1.43,  1.50,  1.41,  1.34,  1.40,
                1.36,  1.30,  1.22,  1.34,  1.27,  1.21,  1.13,  1.15,  1.03,  1.02,  0.98,  0.99,  0.89,  0.85,  0.74,  0.66,  0.50,
                0.41,  0.45,  0.38,  0.26,  0.22,  0.08,  0.10};
        private double[] thresholds = {4, 6, 8, 10, 12, 14};
        private int[][] maxCaContacts = new int[thresholds.length][20];
        int[][] maxCbContacts = new int[thresholds.length][20];
        private MinEstimator[][] minCaDistance = new MinEstimator[20][20];
        private MinEstimator[][] minCbDistance = new MinEstimator[20][20];
        private int[] torsionHistogram = new int[360];
        private int[] torsionHistogramPro = new int[360];
        private int[] torsionHistogramGly = new int[360];
        private AngleStatistics[][][][][][] phiPsi = new AngleStatistics[36][36][8][2][2][2];
        private JSONObject residues;

        public enum LocalPairs {O_CB, O_C, O_O, O_N, C_CB, C_C, C_O, C_N, N_O, N_N, CB_O, CB_N, CB_CB}

        public enum GPTypes {GLY, PRO, OTHER}

        private MinEstimator[][] minLocal = new MinEstimator[GPTypes.values().length][LocalPairs.values().length];

        public enum Neighbors {CA_CA_min, CA_CA_max, CA_CA_minProline, CA_CA_maxProline}

        private MinEstimator[] neighbors = new MinEstimator[Neighbors.values().length];
        private EnergyStatistics energyStatisticsE = new EnergyStatistics();

        public ContactConstraintsAnalyzer() {
            for (int i = 0; i < 20; i++) {
                for (int j = 0; j < 20; j++) {
                    minCaDistance[i][j] = new MinEstimator();
                    minCbDistance[i][j] = new MinEstimator();
                }
            }
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < minLocal[0].length; j++) {
                    minLocal[i][j] = new MinEstimator();
                }
            }
            for (int iTorsion1 = 0; iTorsion1 < 36; iTorsion1++)
                for (int iTorsion2 = 0; iTorsion2 < 36; iTorsion2++)
                    for (int iAngle2 = 0; iAngle2 < 8; iAngle2++)
                        for (int iAngle1 = 0; iAngle1 < 2; iAngle1++)
                            for (int iAngle3 = 0; iAngle3 < 2; iAngle3++)
                                for (int iPhiPsi = 0; iPhiPsi < 2; iPhiPsi++)
                            phiPsi[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][iPhiPsi] = new AngleStatistics();
        }

        public void addResidues(JSONObject residues) {
            this.residues = residues;
        }
        public void toJason(String fileName) throws IOException {
            double[][][][][][][] phiPsiOut = new double[36][36][8][2][2][2][2];
            for (int iTorsion1 = 0; iTorsion1 < 36; iTorsion1++)
                for (int iTorsion2 = 0; iTorsion2 < 36; iTorsion2++)
                    for (int iAngle2 = 0; iAngle2 < 8; iAngle2++) {
                        double[][] last = {{-9999,-9999},{-9999,-9999}};
                        for (int iAngle1 = 0; iAngle1 < 2; iAngle1++)
                            for (int iAngle3 = 0; iAngle3 < 2; iAngle3++)
                                for (int i = 0; i < 2; i++) {
                                    if (!phiPsi[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i].invalid()) {
                                        last[i][0] = 180 * phiPsi[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i].mean() / Math.PI;
                                        phiPsiOut[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i][0] = last[i][0];
                                        last[i][1] =  phiPsi[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i].std();
                                        phiPsiOut[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i][1] = last[i][1];

                                    } else {
                                        phiPsiOut[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i][0] = last[i][0];
                                        phiPsiOut[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][i][1] = last[i][1];
                                    }
                                }
                        }
            MeshiWriter writer = new MeshiWriter(fileName);
            JSONObject jsonObject = new JSONObject();
            jsonObject.put("ResidueTypes", getResidueTypes());
            jsonObject.put("LocalPairs", getLocalPairs());
            jsonObject.put("GPTypes", getGPTypes());
            jsonObject.put("minCaDistance", getMinMatrix(minCaDistance,100));
            jsonObject.put("minCbDistance", getMinMatrix(minCbDistance,100));
            jsonObject.put("minLocal", getMinMatrix(minLocal,100));
            jsonObject.put("alphaHistogram", torsionHistogram);
            jsonObject.put("alphaHistogramGly", torsionHistogramGly);
            jsonObject.put("alphaHistogramPro", torsionHistogramPro);
            jsonObject.put("IndicatorEL", energyStatisticsE.toJason());
            jsonObject.put("phiPsi",phiPsiOut);
            jsonObject.put("residues",residues);
            writer.print(jsonObject);
            writer.close();
        }



        public double[][] getMinMatrix(MinEstimator[][] minEstimators, int threshold) {
            int n = minEstimators.length;
            int m = minEstimators[0].length;
            double[][] out = new double[n][m];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++)
                    out[i][j] = minEstimators[i][j].min(threshold);
            }
            return out;
        }


        public JSONObject getGPTypes() {
            JSONObject out = new JSONObject();
            GPTypes[] types = GPTypes.values();
            for (GPTypes type : types)
                out.put(type.toString(), type.ordinal());
            return out;
        }

        public JSONObject getLocalPairs() {
            JSONObject out = new JSONObject();
            LocalPairs[] pairs = LocalPairs.values();
            for (LocalPairs pair : pairs)
                out.put(pair.toString(), pair.ordinal());
            return out;
        }

        public JSONObject getResidueTypes() {
            JSONObject out = new JSONObject();
            ResidueType[] types = ResidueType.values();
            for (int i = 0; i < 20; i++)
                out.put(types[i].nameOneLetter(), i);
            return out;
        }

        public static double[][] newTable(int n, int m) {
            double[][] out = new double[n][m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    out[i][j] = 1000;
            return out;
        }

        public ProteinInfo analyze(Protein protein) {
            ArrayList<MeshiInfo> infoList = new ArrayList();
            ProteinInfo info = null;
            EnergyInfoElement infoElement;

            Chain chain = protein.chain();
            for (Atom atom : chain.atoms())
                if ((!atom.nowhere()) && (atom.x() < -900))
                    atom.setNowhere();
            int chainLength = chain.size();
            Residue prevResidue = null;
            MolecularSystem ms = chain.firstNonDummyResidue().ca().molecularSystem;
            ms.terminator().reset();
            DistanceMatrix distanceMatrix = new DistanceMatrix(ms);

            for (int i = 1; i < chainLength; i++) {
                analyzeTorsions(chain, i, distanceMatrix);
                Residue residueI = chain.residueAt(i);
                if (!residueI.dummy() && (!residueI.ca().nowhere())) {
                    Atom caI = residueI.ca();
                    Atom cbI = residueI.cb();
                    if (cbI == null) {
                        if (residueI.type == ResidueType.GLY)
                            cbI = caI;
                        else throw new RuntimeException("This is weird.");
                    }
                    ResidueType typeI = residueI.type();
                    if (prevResidue != null)
                        getMinLocal(residueI, prevResidue);
                    for (int j = i + 1; j < chainLength; j++) {
                        Residue residueJ = chain.residueAt(j);
                        if ((!residueJ.dummy()) && (!residueJ.nowhere())) {
                            Atom caJ = residueJ.ca();
                            Atom cbJ = residueJ.cb();
                            if (cbJ == null) {
                                if (residueJ.type == ResidueType.GLY)
                                    cbJ = caJ;
                                else throw new RuntimeException("This is weird.");
                            }
                            ResidueType typeJ = residueJ.type;
                            int separation = residueJ.number() - residueI.number();
                            double distanceCa = (new FreeDistance(caI, caJ)).distance();
                            double distanceCb;
                            if ((!cbI.nowhere()) & (!cbJ.nowhere())) {
                                distanceCb = (new FreeDistance(cbI, cbJ)).distance();
                                if (distanceCb == 0) {
                                    distanceCb = 1000; //throw new RuntimeException("Weird distance in "+protein.name()+"\n"+cbI+"\n"+cbJ);
                                }
                            } else distanceCb = 1000;
                            if (separation == 1) {
                                //getNeighborsData(distanceCa, typeI, typeJ);
                            } else {
                                getData(distanceCa, typeI, typeJ, minCaDistance);
                                getData(distanceCb, typeI, typeJ, minCbDistance);
                            }
                        }
                    }
                    prevResidue = residueI;
                }
            }
            return info;
        }

        private void getMinLocal(Residue residue, Residue prevResidue) {
            Atom ca1 = prevResidue.ca();
            Atom cb1 = prevResidue.ca();
            Atom o1 = prevResidue.carbonylO();
            Atom c1 = prevResidue.carbonylC();
            Atom n1 = residue.amideN();
            Atom ca2 = residue.ca();
            Atom cb2 = residue.cb();
            Atom c2 = residue.carbonylC();
            Atom o2 = residue.carbonylO();
            Atom n2 = residue.amideN();
            int type;
            if (residue.type == ResidueType.GLY)
                type = GPTypes.GLY.ordinal();
            else if (residue.type == ResidueType.PRO)
                type = GPTypes.PRO.ordinal();
            else
                type = GPTypes.OTHER.ordinal();

            double distance;
            if ((!ca1.nowhere()) & (!ca2.nowhere()) & (ca1.distanceFrom(ca2) > 3.75)) {
                if (!n1.nowhere() & (!o2.nowhere())) {
                    distance = n1.distanceFrom(o2);
                    minLocal[type][LocalPairs.N_O.ordinal()].add(distance);
                }
                if ((!n1.nowhere()) & (!n2.nowhere())) {
                    distance = n1.distanceFrom(n2);
                    minLocal[type][LocalPairs.N_N.ordinal()].add(distance);
                }
                if (!o1.nowhere() & (!c2.nowhere())) {
                    distance = o1.distanceFrom(c2);
                    minLocal[type][LocalPairs.O_C.ordinal()].add(distance);
                }
                if ((cb1 != null) && (!cb1.nowhere()) & ((cb2 != null) && (!cb2.nowhere()))) {
                    distance = cb1.distanceFrom(cb2);
                    minLocal[type][LocalPairs.CB_CB.ordinal()].add(distance);
                }
                if ((cb2 != null) && (!cb2.nowhere()) & (!o1.nowhere())) {
                    distance = o1.distanceFrom(cb2);
                    minLocal[type][LocalPairs.O_CB.ordinal()].add(distance);
                }

                if ((!o1.nowhere()) & (!n2.nowhere())) {
                    distance = o1.distanceFrom(n2);
                    minLocal[type][LocalPairs.O_N.ordinal()].add(distance);
                }

                if ((!o1.nowhere()) & (!o2.nowhere())) {
                    distance = o1.distanceFrom(o2);
                    minLocal[type][LocalPairs.O_O.ordinal()].add(distance);
                }

                if ((!c1.nowhere() & (!c2.nowhere()))) {
                    distance = c1.distanceFrom(c2);
                    minLocal[type][LocalPairs.C_C.ordinal()].add(distance);
                }

                if ((cb2 != null) && (!cb2.nowhere()) & (!c1.nowhere())) {
                    distance = c1.distanceFrom(cb2);
                    minLocal[type][LocalPairs.C_CB.ordinal()].add(distance);
                } else minLocal[type][LocalPairs.C_CB.ordinal()].add(1000);

                if ((c1 != null) && (!c1.nowhere()) & (!n2.nowhere())) {
                    distance = c1.distanceFrom(n2);
                    minLocal[type][LocalPairs.C_N.ordinal()].add(distance);
                }

                if ((c1 != null) && (!c1.nowhere()) & (!o2.nowhere())) {
                    distance = c1.distanceFrom(o2);
                    minLocal[type][LocalPairs.C_O.ordinal()].add(distance);
                }

                if ((cb2 != null) && (!cb2.nowhere()) & (!n2.nowhere())) {
                    distance = cb2.distanceFrom(n2);
                    minLocal[type][LocalPairs.CB_N.ordinal()].add(distance);
                } else minLocal[type][LocalPairs.CB_N.ordinal()].add(1000);

                if ((cb2 != null) && (!cb2.nowhere()) & (!o2.nowhere())) {
                    distance = cb2.distanceFrom(o2);
                    minLocal[type][LocalPairs.CB_O.ordinal()].add(distance);
                } else minLocal[type][LocalPairs.CB_O.ordinal()].add(1000);
            }
        }

        private void getData(double distance, ResidueType typeI, ResidueType typeJ, MinEstimator[][] table) {
            table[typeI.ordinal()][typeJ.ordinal()].add(distance);
        }

        private static class MinEstimator {
            private static double MAX_DISANCE = 10;
            private static double BIN_SIZE = 0.1;
            double[] hist;

            public MinEstimator() {
                hist = new double[(int) Math.ceil(MAX_DISANCE / BIN_SIZE)];
            }

            public void add(double distance) {
                if (distance < MAX_DISANCE) {
                    int index = (int) Math.floor(distance / BIN_SIZE);
                        hist[index]++;
                }
            }

            public double min(int threshold) {
                for (int i = 0; i < hist.length; i++) {
                    if (hist[i] >= threshold)
                        return i * BIN_SIZE;
                }
                return 1000;
            }
        }

        public void analyzeTorsions(Chain chain, int i, DistanceMatrix distanceMatrix) {
            double alpha = -9999, prevAlpha = -9999;
            if ((i >= 3) & ( i < chain.size()-2)){
                double phi, psi;
                Residue residueIm2 = chain.residueAt(i-2);
                Residue residueIm1 = chain.residueAt(i-1);
                Residue residueI   = chain.residueAt(i);
                Residue residueIp1 = chain.residueAt(i+1);
                Residue residueIp2 = chain.residueAt(i+2);
                if ((!residueIm2    .dummy()) && (!residueIm1.dummy()) && (!residueI.dummy()) && (!residueIp1.dummy())&&(!residueIp2.dummy())){
                    Atom caIm2 = residueIm2.ca();
                    Atom caIm1 = residueIm1.ca();
                    Atom cIm1 = residueIm1.carbonylC();
                    Atom caI   = residueI.ca();
                    Atom nI    = residueI.amideN();
                    Atom cbI    = residueI.cb();
                    Atom cI    = residueI.carbonylC();
                    Atom caIp1 = residueIp1.ca();
                    Atom nIp1 = residueIp1.amideN();
                    Atom caIp2   = residueIp2.ca();
                    if ((!cI.nowhere()) & (!cIm1.nowhere()) & (!nI.nowhere()) & (!nIp1.nowhere())) {
                        Angle phiAngl1 = new Angle(cIm1, nI, caI, distanceMatrix);
                        Angle phiAngl2 = new Angle(nI, caI, cI, distanceMatrix);
                        phi = (new Torsion(phiAngl1, phiAngl2, distanceMatrix)).torsion();
                        Angle psiAngl1 = new Angle(nI, caI, cI, distanceMatrix);
                        Angle psiAngl2 = new Angle(caI, cI, nIp1, distanceMatrix);
                        psi = (new Torsion(psiAngl1, psiAngl2, distanceMatrix)).torsion();
                        Angle angle0 = new Angle(caIm2, caIm1, caI, distanceMatrix);
                        Angle angle1 = new Angle(caIm1, caI, caIp1, distanceMatrix);
                        Angle angle2 = new Angle(caI, caIp1, caIp2, distanceMatrix);
                        Angle nAngle1 = new Angle(caIm1, caI, nI, distanceMatrix);
                        Angle nAngle2 = new Angle(caIp1, caI, nI, distanceMatrix);
                        Angle cAngle1 = new Angle(caIm1, caI, cI, distanceMatrix);
                        Angle cAngle2 = new Angle(caIp1, caI, cI, distanceMatrix);
                        Angle cbAngle1 = null, cbAngle2 = null;
                        Torsion cbTorsion = null;
                        if (cbI != null) {
                            cbAngle1 = new Angle(caIm1, caI, cbI, distanceMatrix);
                            cbAngle2 = new Angle(caIp1, caI, cbI, distanceMatrix);
                            cbTorsion = new Torsion(cbAngle1, cbAngle2, distanceMatrix);
                        }
                        prevAlpha = (new Torsion(angle0, angle1, distanceMatrix)).torsion();
                        alpha = (new Torsion(angle1, angle2, distanceMatrix)).torsion();
                        int iTorsion1 = (int) Math.floor(18 * (prevAlpha + Math.PI) / Math.PI);
                        int iTorsion2 = (int) Math.floor(18 * (alpha + Math.PI) / Math.PI);
                        int iAngle2 = (int) Math.floor(18 * angle1.angle() / Math.PI) - 8;
                        if (iAngle2 < 0) {
                            iAngle2 = 0;
                            //Utils.printDebug(this, "xxxx " + 180 * angle1.angle() / Math.PI);
                        }
                        if (iAngle2 > 7) {
                            iAngle2 = 7;
                         }
                        int iAngle1 = (int) Math.floor(2 * angle0.angle() / Math.PI);
                        if (iAngle1 > 1) iAngle1 = 1;
                        int iAngle3 = (int) Math.floor(2 * angle2.angle() / Math.PI);
                        if (iAngle3 > 1) iAngle3 = 1;
                            if (residueI.type == ResidueType.GLY)
                                torsionHistogramGly[iTorsion2]++;
                            else if (residueI.type == ResidueType.PRO)
                                torsionHistogramPro[iTorsion2]++;
                            else {
                                torsionHistogram[iTorsion2]++;
                                energyStatisticsE.add(alpha);
                                phiPsi[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][0].add(phi);
                                phiPsi[iTorsion1][iTorsion2][iAngle2][iAngle1][iAngle3][1].add(psi);
                            }
                    }
                }
            }
        }
    }

    private static class IndicatorStatistics extends Statistics{
        double max, min;
        int indicators, all;
        double bottom, top;
        public IndicatorStatistics(double bottom, double top) {
            this.bottom = bottom;
            this .top = top;
            max = indicators = all = n = 0;
            min = 1000;
        }
        public void add(double torsion) {
            double torsionDeg = 180*torsion/Math.PI;
            if ((torsionDeg < top) & (torsionDeg > bottom))
                indicators++;
            all++;
        }

        public void doneProtein(Protein protein) {
            double fraction = (1.0*indicators)/all;
            sum += fraction;
            sum2 += fraction*fraction;
            indicators = all = 0;
            n++;
            if (fraction > max)
                max = fraction;
            if (fraction < min)
                min = fraction;
        }
        public double max() {return max;}
        public double min() {return min;}

        public JSONObject toJason() {
            JSONObject out = super.toJason();
            out.put("max", max());
            out.put("min", min());
            return out;
        }
    }
    private static class EnergyStatistics extends Statistics{
        double max, min;
        double energy;
        int all;
        public EnergyStatistics() {
            energy = all = n = 0;
            max = -1000;
            min = 1000;
        }
        public void add(double torsion) {
            double torsionDeg = 180*torsion/Math.PI;
            int index = (int) Math.floor(180+torsionDeg);
            energy += ENERGY[index];
            all++;
        }

        public void doneProtein(Protein protein) {
            if (energy > 0)
                Utils.printDebug(this,"--------------------------- "+protein+" "+energy);
            energy = energy/all;
            super.add(energy);
            if (energy > max)
                max = energy;
            if (energy < min)
                min = energy;
            energy = 0;
            all = 0;
        }
        public double max() {return max;}
        public double min() {return min;}

        public JSONObject toJason() {
            JSONObject out = super.toJason();
            out.put("max", max());
            out.put("min", min());
            return out;
        }
    }
    private static class Statistics {
        double sum, sum2;
        int n;
        public Statistics() {
            sum = sum2 = n = 0;
        }
        public void add(double torsion) {
            sum += torsion;
            sum2 += torsion*torsion;
            n++;
        }

        public double mean() {
            return sum/n;
        }
        public double std() {
            return Math.sqrt(sum2/n - mean()*mean());
        }
        public JSONObject toJason() {
            JSONObject out = new JSONObject();
            out.put("mean",mean());
            out.put("std",std());
            out.put("n",n);
            return out;
        }
    }
    private static class AngleStatistics extends Statistics {
        private double sumSin, sumCos;
        public AngleStatistics() {
            sumSin = sumCos = 0;
        }
        public void add(double angle) {
            sumSin += Math.sin(angle);
            sumCos += Math.cos(angle);
            n++;
        }
        public double mean() {
            return Math.atan2(sumSin/n, sumCos/n);
        }
        public int n() {return n;}
        public double r2() {
            double meanSin = sumSin/n;
            double meanCos = sumCos/n;
            double r2 = meanSin*meanSin + meanCos*meanCos;
            if (r2 > 1) //numerical noise
                r2 = 1;
            return r2;
        }
        public double std() {
            if (n == 0) return -9999;
            double r2 = r2();
            double std = Math.sqrt(Math.log(1/r2));
            if (Double.isNaN(std))
                throw new RuntimeException("This is weird");
            return std;
        }

        public boolean invalid() {
            return (n() == 0) | (r2() < 0.0001);
        }

    }

    public static JSONObject getResidues(Protein protein) {
        JSONObject residues = new JSONObject();
        ResidueType type = ResidueType.ALA;
        String[] atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.CYS;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.ASP;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.GLU;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.HIS;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.ILE;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.LYS;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.LEU;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.MET;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.ASN;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.PRO;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.GLN;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.ARG;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.SER;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.THR;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.VAL;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.TRP;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        type = ResidueType.TYR;
        atm = getResidue(protein, type);
        if (atm == null)
            throw new RuntimeException("This is weird");
        residues.put(type.toString(),atm);

        return residues;
    }

    private static String[] getResidue(Protein protein, ResidueType type) {
        String[] out = null;
        for (Residue residue : protein.residues())
            if (residue.type == type) {
                out = new String[residue.getAtoms().size()];
                for (int i = 0; i < residue.getAtoms().size(); i++)
                    out[i] = residue.getAtoms().get(i).toString();
                return out;
            }
        return null;
    }

    private static void addSidueChains(Protein protein, JSONObject data) {
        MolecularSystem molecularSystem = new MolecularSystem();
    }
}



