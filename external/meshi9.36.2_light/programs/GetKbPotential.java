package programs;

import meshi.energy.KB2013.*;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 27/07/13
 * Time: 19:22
 * To change this template use File | Settings | File Templates.
 */
public class GetKbPotential{
    public static void main(String argv[]) throws IOException,ParserConfigurationException,SAXException{
        Separator separator = new SeparatorBasic();
        KbStatistics kbs = new KbStatistics(separator, 20, 0.1);
        kbs.read(argv[0]);
        KbPotential kbPotential = new KbPotential(kbs);
        AtomList atoms   = new AtomList(argv[1]);
        Protein protein = new Protein(atoms, new ResidueExtendedAtomsCreator());
        System.out.println(protein);
        kbPotential.add(protein);
//        kbPotential.save(argv[0]+".potential");
    }
}
