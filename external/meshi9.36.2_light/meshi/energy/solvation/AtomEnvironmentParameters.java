package meshi.energy.solvation;

import meshi.energy.Parameters;
import meshi.parameters.AtomType;
import meshi.util.file.MeshiWriter;
import org.w3c.dom.Document;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 19/02/2016.
 */
public class AtomEnvironmentParameters implements Parameters {
    private ArrayList<AtomTypeEnvironmentParameters> types;
    public AtomEnvironmentParameters(File xmlFile, String functionName) {
        types = new ArrayList();
        Document parametersData = readXmlFile(xmlFile);
        NodeList nodes = parametersData.getElementsByTagName("atomEnvironment");
        if (nodes.getLength() == 0)
            throw new RuntimeException("failed to read file "+xmlFile);
        Node typeNode = null;
        for (int iNode = 0; iNode < nodes.getLength(); iNode++) {
            typeNode = nodes.item(iNode);
            if (typeNode != null) {
                types.add(new AtomTypeEnvironmentParameters(typeNode, functionName));
            }
        }
    }

    public AtomTypeEnvironmentParameters getParameters(AtomType atomType) {
        for (AtomTypeEnvironmentParameters atomTypeEnvironmentParameters : types) {
            if (atomTypeEnvironmentParameters.type == atomType) return atomTypeEnvironmentParameters;
        }
        return null;
    }


    private Document readXmlFile(File xmlFile) {
        DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        DocumentBuilder dBuilder = null;
        try {
            dBuilder = dbFactory.newDocumentBuilder();
        } catch (ParserConfigurationException e) {
            System.out.println("**********************************************\nProblem # 1 in ConfigurationArray constructor.\n" + e.getMessage());
            e.printStackTrace();
            throw new RuntimeException("Quitting ConfigurationArray constructor");
        }
        Document out = null;
        try {
            out = dBuilder.parse(xmlFile);
        } catch (SAXException e) {
            System.out.println("**********************************************\nProblem #2 in ConfigurationArray constructor.\n" + e.getMessage());
            e.printStackTrace();
            throw new RuntimeException("Quitting ConfigurationArray constructor");
        } catch (IOException e) {
            System.out.println("**********************************************\nProblem #3 in ConfigurationArray constructor.\n" + e.getMessage());
            e.printStackTrace();
            throw new RuntimeException("Quitting ConfigurationArray constructor");
        }
        out.getDocumentElement().normalize();
        return out;
    }


    public static void main(String[] argv) throws IOException{
        AtomEnvironmentParameters atomEnvironmentParameters = new AtomEnvironmentParameters(new File(argv[0]), "energy");
        AtomTypeEnvironmentParameters matrix = atomEnvironmentParameters.getParameters(AtomType.AN);
        MeshiWriter writerX = new MeshiWriter("testX.txt");
        MeshiWriter writerY = new MeshiWriter("testY.txt");
        MeshiWriter writerZ = new MeshiWriter("testZ.txt");
        for (double cnc = 0; cnc < 20; cnc += 0.1) {
            for (double hbc = 0; hbc < 10; hbc += 0.1) {
                matrix.calc1(cnc, hbc);
                writerX.print(cnc + "\t");
                writerY.print(hbc + "\t");
                writerZ.print(matrix.getValue() + "\t");
            }
            writerX.println();
            writerY.println();
            writerZ.println();
        }
        writerX.close();
        writerY.close();
        writerZ.close();

        matrix.test(0, 0.7);
        System.out.println("-------------");
        matrix.test(0.5, 0);
        System.out.println("-------------");
        matrix.test(0, 0.7);
        System.out.println("-------------");
        matrix.test(0, 0);
        System.out.println("-------------");
        matrix.test(0, 0.5);
        System.out.println("-------------");
        matrix.test(0, 0.25);
        System.out.println("-------------");
        matrix.test(0.5, 0);
        System.out.println("-------------");
        matrix.test(0, 0.4);
        System.out.println("-------------");
        matrix.test(0.5, 0.4);
    }


}
