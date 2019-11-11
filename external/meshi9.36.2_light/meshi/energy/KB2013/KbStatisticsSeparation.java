package meshi.energy.KB2013;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;
import meshi.util.UpdateableException;
import meshi.util.info.MeshiInfoXMLwriter;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;

/**
 *
 */
public class KbStatisticsSeparation {
    public final double maxDistance;
    public final double binSize;
    private KbMatrix matrix;
    private KbVector referenceVector;
    public final int ID;
    public KbStatisticsSeparation(double maxDistance, double binSize,int ID) {
        this.maxDistance = maxDistance;
        this.binSize     = binSize;
        matrix           = new KbMatrix(maxDistance,binSize);
        referenceVector  = new KbVector(maxDistance,binSize);
        this.ID          = ID;
    }


    public void read(org.w3c.dom.Document doc) throws IOException,ParserConfigurationException,SAXException{
        NodeList elements   = doc.getElementsByTagName("KBstatisticsSeparation");
        Node statisticsNode = null;
        for (int iElement = 0; iElement <elements.getLength(); iElement++) {
            Node element = elements.item(iElement);
            if (element.getAttributes().getNamedItem("ID").getTextContent().equals(""+ID))
                statisticsNode = element;
        }
        if (statisticsNode == null)
            throw new RuntimeException("This is weird");
        NodeList children = statisticsNode.getChildNodes();
        Node referenceNode = pickChildByTag(children, "reference");
        fillReferenceVector(referenceNode);
        Node matrixNode = pickChildByTag(children, "KBmatrix");
        matrix.readXML(matrixNode);
    }
    private void fillReferenceVector (Node referenceNode){
        NodeList elements = referenceNode.getChildNodes();
        Node vectorNode = pickChildByTag(elements,"KBvector");
        if (vectorNode == null ) throw new RuntimeException("Cannot fined "+"KBvector");
        referenceVector.readXML(vectorNode);
    }

    protected static Node pickChildByTag(NodeList elements,String tag){
        Node node = null;
        for (int i = 0; i < elements.getLength(); i++) {
            Node nodeI = elements.item(i);
            if (nodeI.getNodeName().equals(tag))
                node = nodeI;
        }
        return node;
    }




    public KbVector getPMF(Atom atomI,Atom atomJ,double weight) {
        KbVector observed = matrix.getData(atomI,atomJ);
        return KbVector.logOds(observed,referenceVector,weight);
    }

    public void writeXml(MeshiInfoXMLwriter writer) throws IOException{
        writer.printf("\t<KBstatisticsSeparation ID=\""+ID+"\">\n");
        writer.printf("\t\t" + "<reference>\n");
        referenceVector.writeXML(writer, "\t\t\t");
        writer.printf("\t\t" + "</reference>\n");
        matrix.writeXML(writer);
        writer.printf("\t</KBstatisticsSeparation>\n");
    }
    protected void add(AtomType type1,AtomType type2, double distance,double weight) throws UpdateableException{
            matrix.addData(type1,type2,distance,weight);
            referenceVector.add(distance, weight);
    }

}
