package meshi.energy.KB2013;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;
import meshi.util.info.MeshiInfoXMLwriter;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

/**
* Created by IntelliJ IDEA.
* User: chen
* Date: 11/07/13
* Time: 17:48
* To change this template use File | Settings | File Templates.
*/
public class KbMatrix {
    private AtomType[] atomTypes = AtomType.values();
    private int nTypes = atomTypes.length;
    private KbVector[][] matrix = new KbVector[nTypes][nTypes];
    private double maxDistance, binSize;
    public KbMatrix(double maxDistance, double binSize) {
        this.maxDistance = maxDistance;
        this.binSize     = binSize;
    }
    public void addData(AtomType type1, AtomType type2,
                        double distance, double weight){
        if (type1.ordinal()<type2.ordinal())
            addData(type2,type1,distance,weight);
        else {
            if (matrix[type1.ordinal()][type2.ordinal()] == null)
                matrix[type1.ordinal()][type2.ordinal()] = new KbVector(maxDistance,binSize);
            matrix[type1.ordinal()][type2.ordinal()].add(distance,weight);
        }
    }
    public KbVector getData(Atom atomI, Atom atomJ){
        AtomType typeI  = atomI.type();
        AtomType typeJ  = atomJ.type();
        if (typeJ.ordinal()>typeI.ordinal())
            return getData(atomJ, atomI);
        else {
            return matrix[typeI.ordinal()][typeJ.ordinal()];
        }
    }

    public void setData(Atom atomI, Atom atomJ, KbVector vector){
        AtomType typeI  = atomI.type();
        AtomType typeJ  = atomJ.type();
        if (typeJ.ordinal()>typeI.ordinal())
            setData(atomJ, atomI,vector);
        else {
            matrix[typeI.ordinal()][typeJ.ordinal()] = vector;
        }
    }

    public void writeXML(MeshiInfoXMLwriter writer) {
        writer.printf("\t\t"+"<KBmatrix>\n");
        AtomType[] types = AtomType.values();
        for (int iType = 0; iType < types.length; iType++)
            for (int jType = 0; jType <= iType; jType++) {
                writer.printf("\t\t\t"+"<KBmatrixEntry iType=\"%s\" jType=\"%s\">\n",
                              types[iType],types[jType]);
                if(matrix[iType][jType] != null)
                    matrix[iType][jType].writeXML(writer,"\t\t\t");
                writer.printf("\t\t\t"+"</KBmatrixEntry>\n");
            }
        writer.printf("\t\t"+"</KBmatrix>\n");
    }

    public void readXML(Node node){
        NodeList elements = node.getChildNodes();
        for (int i = 0; i< elements.getLength(); i++) {
            Node item = elements.item(i);
            if (item.getNodeName().equals("KBmatrixEntry")) {
                NamedNodeMap attributes = item.getAttributes();
                Node iType = attributes.getNamedItem("iType");
                Node jType = attributes.getNamedItem("jType");
                AtomType typeI = AtomType.type(iType.getNodeValue());
                AtomType typeJ = AtomType.type(jType.getNodeValue());
                NodeList itemElements = item.getChildNodes();
                Node vectorNode = KbStatisticsSeparation.pickChildByTag(itemElements, "KBvector");
                if (vectorNode != null) {
                    matrix[typeI.ordinal()][typeJ.ordinal()] = new KbVector(maxDistance,binSize);
                    KbVector vector = matrix[typeI.ordinal()][typeJ.ordinal()];
                    vector.readXML(vectorNode);
                }
            }
         }
    }

}
