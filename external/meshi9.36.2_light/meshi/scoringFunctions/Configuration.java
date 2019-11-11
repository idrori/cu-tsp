package meshi.scoringFunctions;

import org.w3c.dom.Node;

/**
* Created with IntelliJ IDEA.
* User: chen
* Date: 18/12/13
* Time: 10:59
* To change this template use File | Settings | File Templates.
*/
class Configuration {
    static Node currentNode;
    public final String[] fields, normalizers;
    public final String objectiveFunction;
    public final double[] exponentValues, coefs;
    public final int[][] exponentIndices;
    public final double sigmoidAlpha, sigmoidBeta, sigmoidGamma, sigmoidDelta;
    String sigmoidString;
    String[] exponentIndicesStrings, exponentDimensionsStrings, coefsStrings, exponentValuesStrings;
    public Configuration(Node configurationNode) {
        fields             = getNode(configurationNode.getFirstChild().getNextSibling(), "fields");
        normalizers        = getNode(currentNode.getNextSibling().getNextSibling(),"normalizers");

        coefsStrings       = getNode(currentNode.getNextSibling().getNextSibling(),"coefs");
        coefs = new double[coefsStrings.length];
        int i = 0;
        for (String s:coefsStrings) {
            coefs[i] = Double.parseDouble(s);
            i++;
        }

        exponentIndicesStrings    = getNode(currentNode.getNextSibling().getNextSibling(),"exponentIndices");
        exponentDimensionsStrings = getNode(currentNode,"exponentIndices","dim");
        exponentIndices = new int[Integer.parseInt(exponentDimensionsStrings[0])][Integer.parseInt(exponentDimensionsStrings[1])];
        int k = 0;
        for (int iDim = 0; iDim < exponentIndices.length; iDim++)
            for (int jDim = 0; jDim < exponentIndices[0].length; jDim++) {
                exponentIndices[iDim][jDim] = Integer.parseInt(exponentIndicesStrings[k]);
                k++;
            }

        exponentValuesStrings     = getNode(currentNode.getNextSibling().getNextSibling(),"exponentValues");
        exponentValues            = new double[exponentValuesStrings.length];
        i = 0;
        for (String s:exponentValuesStrings) {
            exponentValues[i] = Double.parseDouble(s);
            i++;
        }



        objectiveFunction  = getNode(currentNode.getNextSibling().getNextSibling(),"objectiveFunction")[0];

        sigmoidString = getNode(currentNode.getNextSibling().getNextSibling(), "sigmoid","alpha")[0];
        sigmoidAlpha  = Double.parseDouble(sigmoidString);
        sigmoidString = getNode(currentNode, "sigmoid","beta")[0];
        sigmoidBeta  = Double.parseDouble(sigmoidString);
        sigmoidString = getNode(currentNode, "sigmoid","gamma")[0];
        sigmoidGamma  = Double.parseDouble(sigmoidString);
        sigmoidString = getNode(currentNode, "sigmoid","delta")[0];
        sigmoidDelta  = Double.parseDouble(sigmoidString);
    }

    private static Node checkNode(Node node,String name) {
        if (!node.getNodeName().equals(name))
            throw new RuntimeException("Weird node "+node.getNodeName()+" "+name);
        currentNode = node;
        return node;
    }

    private String[] getNode(Node node,String name) {
        return getNode(node, name,name);
    }
    private String[] getNode(Node node,String nodeName, String attributeName) {
        node = checkNode(node,nodeName);
        return node.getAttributes().getNamedItem(attributeName).getNodeValue().split("\\s+");
    }

 }
