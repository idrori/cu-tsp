package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPair;
import meshi.molecularElements.atoms.AtomPairList;
import meshi.molecularElements.extendedAtoms.ExtendedAtomsProtein;
import meshi.util.Utils;
import meshi.util.filters.IsCSG;

import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 14/06/12
 * Time: 14:57
 * To change this template use File | Settings | File Templates.
 */
public class DetectDisulfide {
    public static void main(String [] argv) throws IOException {
        Utils.verboseOff();
        double[] thresholds = new double[argv.length];
        int expectedN = argv.length;
        double maxThreshold = -100;
        for (int i = 0; i < expectedN; i++) {
                thresholds[i] = Double.valueOf(argv[i]);
                if (maxThreshold < thresholds[i]) maxThreshold = thresholds[i];
        }
        File dir = new File(".");
        File[] files = dir.listFiles();

        int counter = 0;
        for (File file:files) {
            boolean OK = true;
            if (file.getName().endsWith("pdb")) {
                Protein model = ExtendedAtomsProtein.getExtendedAtomsProteinFromPdbFile(file);
                AtomList atoms = model.atoms();
                AtomPairList disulfides = Utils.detectContacts(atoms, IsCSG.isCSG,maxThreshold);
                if (disulfides.size() == expectedN) {
                    for (int i = 0; (i < expectedN) && OK; i++) {
                        AtomPair disulfide  = disulfides.get(i);
                        if (disulfide.atom1().distanceFrom(disulfide.atom2())> thresholds[i])
                            OK = false;
                    }
                    if (OK) {
                       System.out.print("\n"+(++counter)+" "+file.getName()+" ");
                        for (int i = 0; i < expectedN; i++) {
                            AtomPair disulfide  = disulfides.get(i);
                            Atom atom1 = disulfide.atom1();
                            Atom atom2 = disulfide.atom2();
                            System.out.print(atom1.residueNumber()+" "+atom1.number()+" "+
                                     atom2.residueNumber()+" "+atom2.number()+" "+atom1.distanceFrom(atom2)+", ");
                        }
                       System.out.println();
                    }
                }
                else OK = false;
                if (!OK) System.out.print('.');
            }
        }
    }

}
