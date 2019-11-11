package meshi.energy.solvation;

import meshi.parameters.AtomType;
import meshi.util.Utils;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.util.ArrayList;

/**
 * Created by chen on 20/02/2016.
 */
public class AtomTypeEnvironmentParameters {
    public static final double ALPHA = 100;
    public static final double EPSILON = 1/(ALPHA+1);
    public final double MAX_VALUE;
    AtomType type;
    double[][] environmentMatrix;
    private double value, dValue, dValueDhbc, dValueDcnc, sumWeights, sumDweights;
    private double[][] values = new double[3][3];
    private double[][] dCnc   = new double[3][3];
    private double[][] dHbc   = new double[3][3];
    private double[][] dis2   = new double[3][3];
    private double[][] weight = new double[3][3];
    private double[][] dWeight = new double[3][3];

    public AtomTypeEnvironmentParameters(Node typeNode, String functionName) {
            Node node = typeNode.getAttributes().item(0);
            String name = node.getNodeValue();
            type = AtomType.type(name);
//            System.out.println(type);
            environmentMatrix = getEnvironmentMatrix(typeNode, functionName);
            value = dValueDcnc = dValueDhbc = Double.MAX_VALUE;
            double max;
            if (environmentMatrix == null) {
                max = 0;
            }
            else {
                max = -10000;
                for (int i = 0; i < environmentMatrix.length; i++)
                    for (int j = 0; j < environmentMatrix[0].length; j++)
                        if (environmentMatrix[i][j] > max) max = environmentMatrix[i][j];
            }
        MAX_VALUE = max;
    }

    public void calc1 (double cnc, double hbc) {
        if (environmentMatrix == null) return;
        if ((cnc >= environmentMatrix.length) || (hbc >= environmentMatrix[0].length / 2)) {
            value = MAX_VALUE;
            return;
        }
        value = dValue = dValueDcnc = dValueDhbc = sumWeights = sumDweights = 0;
        double gama     = 100;
        double epsilon  = 1/(gama+1);
        double maxAlpha = environmentMatrix.length - 1;
        double maxBeta  = environmentMatrix[0].length - 1;
        double[] V      = new double[5];
        double[] dVdCnc = new double[5];
        double[] dVdHbc = new double[5];
        double[] W      = new double[5];
        double[] dWdCnc = new double[5];
        double[] dWdHbc = new double[5];
        double[] alpha  = new double[5];
        double[] beta   = new double[5];
        double[] A      = new double[5];
        double lorenz, Q, dQdCnc, dQdHbc;

        for (int i = 0; i < 5; i++) {
            V[i]      = 0;
            dVdCnc[i] = 0;
            dVdHbc[i] = 0;
            W[i]      = 0;
            dWdCnc[i] = 0;
            dWdHbc[i] = 0;
            alpha[i]  = 0;
            beta[i]   = 0;
            A[i]      = 0;
        }
        Q      = 0;
        dQdCnc = 0;
        dQdHbc = 0;

        int iCnc = (int) Math.floor(cnc - 0.5) + 1; // 0 => 0; 0.5 = > 1; 1 => 1
        int jHbc = (int) Math.floor(2 * hbc - 0.5) + 1; // 0 => 0; 0.5  => 1 ; 1.5 => 3
        int iElelment = 0;
        if ((iCnc >= environmentMatrix.length) || (jHbc >= environmentMatrix[0].length)) {
            value += MAX_VALUE;
            dValueDcnc = 0;
            dValueDhbc = 0;
            return;
        }
        int iStart = (iCnc == 0) ? 0 : -1;
        int iEnd   = (iCnc < environmentMatrix.length-1) ? 2 : 1;
        int jStart = (jHbc == 0) ? 0 : -1;
        int jEnd   = (jHbc < environmentMatrix[0].length-1) ? 2 : 1;
        for (int i = iStart; i < iEnd; i++) {
            double da    = cnc - (iCnc + i);
            if (da >= 1) continue;
            for (int j = jStart; j < jEnd; j++) {
                double db    = hbc - (jHbc + j)/2.0;
                if (db >= 0.5) continue;
                double dadb  = da * da + 4 * db * db;
                if (dadb > 1) continue;
                if (iElelment >= 5)
                    throw new RuntimeException("This is weird.");
                V[iElelment]      = gama * (da * da + 4 * db * db) + 1;
                dVdCnc[iElelment] = 2 * gama * da;
                dVdHbc[iElelment] = 2 * gama * 4*db;
                lorenz            = 1 / V[iElelment];
                W[iElelment]      = (lorenz - epsilon) * (lorenz - epsilon);
                dWdCnc[iElelment] = -2 * (lorenz - epsilon) * dVdCnc[iElelment] / (V[iElelment] * V[iElelment]);
                dWdHbc[iElelment] = -2 * (lorenz - epsilon) * dVdHbc[iElelment] / (V[iElelment] * V[iElelment]);
                Q      += W[iElelment];
                dQdCnc += dWdCnc[iElelment];
                dQdHbc += dWdHbc[iElelment];
                A[iElelment] = environmentMatrix[iCnc + i][jHbc + j];
                iElelment++;
            }
        }

        double iElementMax = iElelment;
        value = 0;
        for (iElelment = 0; iElelment < iElementMax; iElelment++) {
            value += A[iElelment] * W[iElelment] / Q;
            dValueDcnc += A[iElelment] * (dWdCnc[iElelment] - W[iElelment]*dQdCnc/Q)/Q;
            dValueDhbc += A[iElelment] * (dWdHbc[iElelment] - W[iElelment]*dQdHbc/Q)/Q;
        }
    }

        public void test(double cnc,double hbc) {
            double d = 0.00001;
            calc1(cnc, hbc);
            double value = getValue();
            System.out.println("cnc & hbc       = " + cnc + " "+hbc);
            System.out.println("Value           = " + value);
            System.out.println("Analitical dCnc = " + getDvalueDcnc());
            System.out.println("Analitical dHbc = " + getDvalueDhbc());
            calc1(cnc + d, hbc);
            double newValue = getValue();
            System.out.println("Numerical dCnc = " + (newValue - value) / d);
            calc1(cnc, hbc + d);
            newValue = getValue();
            System.out.println("Numerical dHbc = " + (newValue - value)/d);
        }


    public double getValue() {
            if (value == Double.MAX_VALUE)
                return 0;
            double out = value;
            value = Double.MAX_VALUE;
            return out;
        }

        public double getDvalueDhbc() {
            if (dValueDhbc == Double.MAX_VALUE)
                throw new RuntimeException("This is weird");
            double out = dValueDhbc;
            dValueDhbc = Double.MAX_VALUE;
            return out;
        }
        public double getDvalueDcnc() {
        if (dValueDcnc == Double.MAX_VALUE)
            throw new RuntimeException("This is weird");
            double out = dValueDcnc;
            dValueDcnc = Double.MAX_VALUE;
        return out;
    }



        private double[][] getEnvironmentMatrix(Node typeNode, String functionName) {
            Node functionNode = getFunctionNode(typeNode, functionName);
            ArrayList<Node> matrixCncNodes = getMatrixRowNodes(functionNode);
            ArrayList<ArrayList<Node>> matrixNodes = getMatrixNodes(matrixCncNodes);
            if (matrixNodes != null) {
                int nCnc = matrixCncNodes.size();
                int nHbc = matrixNodes.get(0).size();
                environmentMatrix = new double[nCnc][nHbc];
                for (int iCnc = 0; iCnc < nCnc; iCnc++)
                    for (int jHbc = 0; jHbc < nHbc; jHbc++) {
                        Node node1 = matrixNodes.get(iCnc).get(jHbc);
                        String valueString = node1.getAttributes().item(1).getNodeValue();
                        environmentMatrix[iCnc][jHbc] = Double.valueOf(valueString);
                    }
            } else {
                environmentMatrix = null;
            }
            return environmentMatrix;
        }

        private ArrayList<ArrayList<Node>> getMatrixNodes(ArrayList<Node> matrixCncNodes) {
            if (matrixCncNodes != null) {
                ArrayList<ArrayList<Node>> matrixNodes = new ArrayList();
                for (Node cncNode : matrixCncNodes) {
                    NodeList hbcNodeList = cncNode.getChildNodes();
                    ArrayList<Node> hbcNodes = new ArrayList();
                    for (int jNode = 0; jNode < hbcNodeList.getLength(); jNode++) {
                        Node hbcNode = hbcNodeList.item(jNode);
                        if (hbcNode.getNodeName().equals("hbc"))
                            hbcNodes.add(hbcNode);
                    }
                    matrixNodes.add(hbcNodes);
                }
                return matrixNodes;
            }
            return null;
        }

        private ArrayList<Node> getMatrixRowNodes(Node functionNode) {
            ArrayList<Node> out = new ArrayList();
            if (functionNode != null) {
                NodeList matrixRowNodes = functionNode.getChildNodes();
                for (int iRowNode = 0; iRowNode < matrixRowNodes.getLength(); iRowNode++) {
                    Node rowNode = matrixRowNodes.item(iRowNode);
                    if (rowNode.getNodeName().equals("cnc")) {
                        out.add(rowNode);
                    }
                }
                return out;
            }
            return null;
        }

        private Node getFunctionNode(Node typeNode, String functionName) {
            Node functionNode = null;
            NodeList functions = typeNode.getChildNodes();
            if (functions.getLength() > 1) {
                for (int iFunction = 0; iFunction < functions.getLength(); iFunction++) {
                    Node node1 = functions.item(iFunction);
                    String name = node1.getNodeName();
                    if (name.equals(functionName)) {
                        functionNode = node1;
                        break;
                    }
                }
                if (functionNode == null)
                    throw new RuntimeException("could not find node for " + functionName);
            }
            return functionNode;
        }
}
