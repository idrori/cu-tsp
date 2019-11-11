package programs;

import meshi.energy.KB2013.AtomPairIterator;
import meshi.energy.KB2013.KbPotential;
import meshi.energy.KB2013.KbStatistics;
import meshi.energy.KB2013.Separator;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPair;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentException;
import meshi.util.Rms;
import meshi.util.file.MeshiWriter;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 28/07/13
 * Time: 07:32
 * To change this template use File | Settings | File Templates.
 */
public class KbTest {
    public static void main(String[] argv) throws IOException,ParserConfigurationException,SAXException,AlignmentException{
        String targetName = argv[0];
        String dbName = argv[1];
        Separator separator = Separator.getSeparator(argv[2]);
        double maxDistance = Double.parseDouble(argv[3]);
        double binSize = Double.parseDouble(argv[4]);
        KbStatistics kbStatistics = new KbStatistics(separator, maxDistance,binSize);
        kbStatistics.read(dbName);
        AtomList atoms   = new AtomList(targetName+"/"+targetName+".N.pdb");
        Protein nativeStructure = new Protein(atoms, new ResidueExtendedAtomsCreator());
        KbPotential potential = new KbPotential(kbStatistics);
        potential.add(nativeStructure);
        File dir = new File(targetName);
        File[] files = dir.listFiles();
        MeshiWriter writer = new MeshiWriter(targetName+"."+dbName.substring(0,2)+".dat");

        int i = 0;
        for (File file : files)  {
            if (file.getName().startsWith(targetName) && (file.getName().indexOf(".N.pdb")<0)) {
                atoms         = new AtomList(file.getAbsolutePath());
                Protein model = new Protein(atoms, new ResidueExtendedAtomsCreator());
                double[] gdt;
                double score;
                try {
                    gdt  = Rms.gdt(nativeStructure,model);
                    model.atoms().molecularSystem.terminator().reset();
                    AtomPairIterator atomPairIterator = new AtomPairIterator(model.atoms(),4);
                    score = 0;
                    while(atomPairIterator.hasNext()) {
                        AtomPair pair = atomPairIterator.next();
                        Atom atom1 = pair.atom1();
                        Atom atom2 = pair.atom2();
                        double distance = atom1.distanceFrom(atom2);
                        if (distance <= maxDistance){
                            score += potential.getScore(atom1,atom2,distance,model);
                        }
                    }
                }
                catch (Exception ex) {
                    System.out.println(ex.getMessage());
                    ex.printStackTrace();
                    continue;
                }
                System.out.printf("%-10.2f  %-10.2f  %-10.2f   %s\n",gdt[0]*100,score,score/model.chain().numberOfNonDummyResidues(),model.name());
                writer.printf("%-10.2f  %-10.2f  %-10.2f\n",gdt[0]*100,score,score/model.chain().numberOfNonDummyResidues());
            }
        }
        writer.close();
    }
}
