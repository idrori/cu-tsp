package meshi.energy.KB2013;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomPair;
import meshi.parameters.AtomType;

/**
 *
 * */
public class KbPotential {
    private KbMatrix[] matrices;
    private KbStatistics statistics = null;
    private Separator separator;

    public KbPotential(KbStatistics statistics) {
        this.separator = statistics.separator;
        this.statistics = statistics;
        matrices = new KbMatrix[separator.numberOfSeparations()];
        for (int i = 0; i < matrices.length; i++)
            matrices[i] = new KbMatrix(statistics.maxDistance,statistics.binSize);
    }


   public void add(Protein protein) {
       AtomPairIterator atomPairIterator = new AtomPairIterator(protein.atoms(),4);

       while(atomPairIterator.hasNext()) {
           AtomPair  pair = atomPairIterator.next();
           Atom     atom1 = pair.atom1();
           Atom     atom2 = pair.atom2();
           int separation = separator.separation(atom1,atom2);
           if (matrices[separation].getData(atom1,atom2)==null)
                matrices[separation].setData(atom1,atom2,statistics.getPMF(atom1,atom2));
       }
    }

    public double getScore(Atom atom1, Atom atom2, double distance, Protein protein) {
        AtomType type1 = atom1.type();
        AtomType type2 = atom2.type();
        if (type1.ordinal()<type2.ordinal())
            return getScore(atom2,atom1,distance, protein);
        int separation = separator.separation(atom1,atom2);
        KbVector vector = matrices[separation].getData(atom1,atom2);
        if (vector == null) {
            if (protein == null)
                throw new RuntimeException("This is weird");
            add(protein);
            System.out.println(protein.name()+" added to the potential");
            return getScore(atom1,atom2,distance,null);
        }
        return vector.getScore(distance);
    }

//    public void save(String fileName) throws IOException{
//        MeshiInfoXMLwriter writer = new MeshiInfoXMLwriter(fileName);
//        writer.printf("<KBpotential>\n");
//        matrix.writeXML(writer);
//        writer.printf("</KBpotential>\n");
//        writer.close();
//    }
//
//    public void read(String fileName) throws IOException,ParserConfigurationException,SAXException {
//        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
//        DocumentBuilder db = dbf.newDocumentBuilder();
//        org.w3c.dom.Document doc = db.parse(new File(fileName));
//        matrix.readXML(doc);
//    }

}