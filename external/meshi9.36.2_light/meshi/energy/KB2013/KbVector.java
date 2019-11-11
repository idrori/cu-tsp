package meshi.energy.KB2013;

import meshi.util.info.MeshiInfoXMLwriter;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

import java.util.StringTokenizer;

/**
* Created by IntelliJ IDEA.
* User: chen
* Date: 11/07/13
* Time: 17:49
* To change this template use File | Settings | File Templates.
*/
public class KbVector {
    public static final double EPSILON = 0.001;
    public final double maxDistance;
    public final double binSize;
    public final double[] data;
    private IndexAndWeights indexAndWeights = new IndexAndWeights();
    public KbVector(double maxDistance, double binSize) {
         this.maxDistance = maxDistance;
         this.binSize     = binSize;
         data = new double[(int)Math.ceil(maxDistance/binSize)+1];
    }

    public static KbVector scalarMultiplication(KbVector inVector,
                                                double multiplier,KbVector outVector) {
        if (inVector.data.length != outVector.data.length)
            throw new RuntimeException("This is weird.");
        for (int i = 0; i < inVector.data.length; i++) {
            outVector.data[i] = inVector.data[i]*multiplier;
        }
        return outVector;
    }

    public static KbVector devideVectors (KbVector numeratorVector, KbVector denominatorVector,
                                          KbVector quotientVector) {
        int vectorLength = numeratorVector.data.length;
        if ((vectorLength != denominatorVector.data.length) ||
            (vectorLength != quotientVector.data.length))
            throw new RuntimeException("This is weird.");
        for (int i = 0; i < vectorLength; i++) {
            quotientVector.data[i] = (numeratorVector.data[i]+EPSILON)/
                                             (denominatorVector.data[i]+EPSILON);
        }
        return quotientVector;
    }
    public static KbVector logOds(KbVector observed, KbVector reference,double weight) {
        if (observed.data.length != reference.data.length)
            throw new RuntimeException("This is weird.");
        KbVector out = new KbVector(observed.maxDistance,observed.binSize);
        out = scalarMultiplication(reference,weight,out);
        out = devideVectors(observed,out,out);
        for (int i = 0; i < out.data.length; i++)
            out.data[i] = -Math.log(out.data[i]);
        return out;
    }

    public void add(double distance,double weight){
        if (distance < binSize/2) distance = binSize/2;
        if (distance >= maxDistance) return;
        indexAndWeights.calc(distance,binSize);
        data[indexAndWeights.index-1] += indexAndWeights.wBelow * weight;
        data[indexAndWeights.index]   += indexAndWeights.wHere  * weight;
        data[indexAndWeights.index+1] += indexAndWeights.wAbove * weight;
    }

    public double getScore(double distance) {
        indexAndWeights.calc(distance,binSize);
        return  data[indexAndWeights.index-1] * indexAndWeights.wBelow +
                data[indexAndWeights.index]   * indexAndWeights.wHere  +
                data[indexAndWeights.index+1] * indexAndWeights.wAbove;
    }

    private static class IndexAndWeights {
        int index;
        double wBelow, wHere, wAbove;
        public void calc(double distance,double binSize) {
            index = (int) Math.floor(distance/binSize);
            if (index < 1) index = 1;
            double binCenter = index*binSize+binSize/2;
            if (distance >= binCenter){
                wBelow = 0;
                wHere  = (binCenter+binSize-distance)/binSize;
                wAbove = (distance-binCenter)/binSize;
            }
            else {
                wBelow = (binCenter-distance)/binSize;
                wHere  = (distance-(binCenter-binSize))/binSize;
                wAbove = 0;
            }
        }
    }

    public void writeXML(MeshiInfoXMLwriter writer,String blanks) {
        writer.printf(blanks+"<KBvector maxDistance=\"%s\" binSize=\"%s\">\n",maxDistance,binSize);
        int i = 0;
        for (double value : data) {
            if (i == 0)
                writer.print(blanks+"          ");
            writer.printf(" %-10.2f",value);
            i++;
            if (i == 10) {
                writer.print("\n");
                i = 0;
            }
        }
        if(i != 0)
            writer.print("\n");
        writer.printf(blanks+"</KBvector>\n");
    }


    public void readXML(Node vectorNode) {
        NamedNodeMap attributes = vectorNode.getAttributes();
        Node maxDistanceNode = attributes.getNamedItem("maxDistance");
        double maxDistance = Double.parseDouble(maxDistanceNode.getNodeValue());
        if (this.maxDistance != maxDistance)
            throw new RuntimeException("maxDistance does not match "+maxDistance+" "+this.maxDistance);
        Node binSizeNode = attributes.getNamedItem("binSize");
        double binSize = Double.parseDouble(binSizeNode.getNodeValue());
        if (this.binSize != binSize)
            throw new RuntimeException("binSize does not match "+binSize+" "+this.binSize);
        String content = vectorNode.getTextContent();
        StringTokenizer tokenizer = new StringTokenizer(content);
        int i = 0;
        while(tokenizer.hasMoreTokens()) {
            String token = tokenizer.nextToken();
            this.data[i] = Double.parseDouble(token);
            i++;
        }
    }

}
