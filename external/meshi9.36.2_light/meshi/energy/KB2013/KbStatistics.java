package meshi.energy.KB2013;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPair;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.AtomType;
import meshi.util.UpdateableException;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.MeshiInfoXMLwriter;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 05/08/13
 * Time: 08:43
 * To change this template use File | Settings | File Templates.
 */
public class KbStatistics  {
    private KbStatisticsSeparation[] kbStatisticsSeparations;
    public final Separator separator;
    private double[] typeCounters;
    private double[] counter;
    public final double maxDistance, binSize;


    public KbStatistics(Separator separator, double maxDistance, double binSize) {
        this.maxDistance = maxDistance;
        this.binSize     = binSize;
        this.separator = separator;
        kbStatisticsSeparations = new KbStatisticsSeparation[separator.numberOfSeparations()];
        for (int i = 0; i < separator.numberOfSeparations();i++)
            kbStatisticsSeparations[i] = new KbStatisticsSeparation(maxDistance,binSize,i);
        typeCounters     = new double[AtomType.values().length];
        counter          = new double[1];
    }

    public void learn(String dataFilesListName) throws IOException,UpdateableException {
        String          dataLine;
        File dataFilesList       = new File(dataFilesListName);
        MeshiLineReader dataFilesListReader = new MeshiLineReader(dataFilesList);
        while ((dataLine = dataFilesListReader.readLine())!= null){
            String pdbFileName = dataLine.substring(0,4).toUpperCase();
            String chainName   = dataLine.substring(4,5);
            double weight      = Double.parseDouble(dataLine.substring(6));
            System.out.printf("%s %s %f\n", pdbFileName, chainName, weight);
            AtomList atoms   = new AtomList("pdb1.5A/"+pdbFileName);
            Protein protein = new Protein(atoms, new ResidueExtendedAtomsCreator());
            Chain chain   = protein.getChain(chainName);
            System.out.println(chain);
            extractStatistics(chain,weight);
        }
    }
    private void extractStatistics(Chain chain, double weight) throws UpdateableException{
        MolecularSystem molecularSystem = chain.atoms().molecularSystem();
        molecularSystem.terminator().reset();
        AtomPairIterator atomPairIterator = new AtomPairIterator(chain.atoms(),4);

        for (Atom atom:chain.atoms())
                if ((!atom.isHydrogen())&& (!atom.nowhere())) {
                    typeCounters[atom.type().ordinal()] += weight;
                    counter[0] += weight;
                }

        while(atomPairIterator.hasNext()) {
            AtomPair pair = atomPairIterator.next();
            Atom atom1 = pair.atom1();
            Atom atom2 = pair.atom2();
            double distance = atom1.distanceFrom(atom2);
            int separationIndex = separator.separation(atom1,atom2);
            kbStatisticsSeparations[separationIndex].add(atom1.type(), atom2.type(),
                    distance, weight);
        }
    }

    public void save(String fileName) throws IOException{
        MeshiInfoXMLwriter writer = new MeshiInfoXMLwriter(fileName);
        writer.printf("<KBstatistics>\n");
        writer.printf("\t<Separator name=\""+separator.name()+"\" />\n");
        writer.printf("\t" + "<totalCount value=\"%f\"/>\n", counter[0]);
        writer.printf("\t" + "<typeCounters>\n");
        for (AtomType atomType :AtomType.values())
            if(typeCounters[atomType.ordinal()] > 0)
                writer.printf("\t\t"+"<%-4s value=\"%-10.2f\"/>\n",
                              atomType.toString(),typeCounters[atomType.ordinal()]);
        writer.printf("\t"+"</typeCounters>\n");
        for (KbStatisticsSeparation kbStatisticsSeparation : kbStatisticsSeparations)
            kbStatisticsSeparation.writeXml(writer);
        writer.printf("</KBstatistics>\n");
        writer.close();
    }

    public void read(String fileName) throws IOException,ParserConfigurationException,SAXException {
        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
        DocumentBuilder db = dbf.newDocumentBuilder();
        org.w3c.dom.Document doc = db.parse(new File(fileName));
        NodeList separatorNodeList = doc.getElementsByTagName("Separator");
        Node     separatorNode = separatorNodeList.item(0);
        NamedNodeMap attributes = separatorNode.getAttributes();
        if (!attributes.getNamedItem("name").getNodeValue().equals(separator.name()))
            throw new RuntimeException("This is weird\n"+
                    attributes.getNamedItem("name").getNodeValue()+"\n"+
                    separator.name());
        fillTypeCounters(doc);
        fillCounter(doc);
        for (KbStatisticsSeparation kbStatisticsSeparation : kbStatisticsSeparations) {
            kbStatisticsSeparation.read(doc);
        }
    }

    private void fillTypeCounters(org.w3c.dom.Document doc) {
        for (AtomType type: AtomType.values()) {
            NodeList element = doc.getElementsByTagName(type.toString());
            if (element.getLength() > 0) {
                NamedNodeMap values = element.item(0).getAttributes();
                Node value = values.getNamedItem("value");
                String  countString = value.getNodeValue();
                typeCounters[type.ordinal()] = Double.parseDouble(countString);
            }
        }
    }

    private void fillCounter(org.w3c.dom.Document doc) {
        NodeList element = doc.getElementsByTagName("totalCount");
        NamedNodeMap values = element.item(0).getAttributes();
        Node          value = values.getNamedItem("value");
        String  countString = value.getNodeValue();
        counter[0]  = Double.parseDouble(countString);
    }

    private double getWeight(Atom atom1, Atom atom2){
        double weight = (typeCounters[atom1.type().ordinal()]/counter[0])*
                 (typeCounters[atom2.type().ordinal()]/counter[0]);
        return weight;
    }

    public KbVector getPMF(Atom atom1,Atom atom2) {
        double weight = getWeight(atom1,atom2);
        int separationIndex = separator.separation(atom1,atom2);
        return kbStatisticsSeparations[separationIndex].getPMF(atom1,atom2,weight);
    }

}
