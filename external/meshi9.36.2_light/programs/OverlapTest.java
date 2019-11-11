package programs;

import meshi.geometry.Coordinates;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.sequences.AtomAlignment;
import meshi.util.MeshiException;
import meshi.util.Utils;
import meshi.util.overlap.Overlap;

import java.io.File;

/**
 *
 */
public class OverlapTest {
    public static void main(String[] argv) throws Exception {
        Coordinates centerOfMass0 = new Coordinates();
        Coordinates centerOfMass1 = new Coordinates();
        double[][] rotateMatrix;
        double rms;
        double[][] coordinates0;
        double[][] coordinates1;
        double centerOfMassX1 = 0, centerOfMassY1 = 0, centerOfMassZ1 = 0;
        double centerOfMassX0 = 0, centerOfMassY0 = 0, centerOfMassZ0 = 0;
        int size0, size1, size;
        String name0, name1;
        AtomList atomList0, atomList1;
        Atom atom0, atom1;

        Protein protein0 = Protein.getCAproteinFromApdbFile(new File(argv[0]));
        Protein protein1 = Protein.getCAproteinFromApdbFile(new File(argv[1]));

        name0 = protein0.name();
        name1 = protein0.name();

        atomList0 = protein0.atoms();
        atomList1 = protein1.atoms();

        size0 = atomList0.size();
        size1 = atomList1.size();
        if (size0 != size1) throw new Exception("Protein sizes must me equal. size0 = " + size0 + " size1 = " + size1);
        size = size0;

        coordinates0 = new double[3][size];
        coordinates1 = new double[3][size];

        for (int i = 0; i < size; i++) {
            atom0 = atomList0.atomAt(i);
            atom1 = atomList1.atomAt(i);
            coordinates0[0][i] = atom0.x();
            coordinates0[1][i] = atom0.y();
            coordinates0[2][i] = atom0.z();
            coordinates1[0][i] = atom1.x();
            coordinates1[1][i] = atom1.y();
            coordinates1[2][i] = atom1.z();

            centerOfMassX0 += atom0.x();
            centerOfMassY0 += atom0.y();
            centerOfMassZ0 += atom0.z();
            centerOfMassX1 += atom1.x();
            centerOfMassY1 += atom1.y();
            centerOfMassZ1 += atom1.z();
        }

        centerOfMassX0 /= size;
        centerOfMassY0 /= size;
        centerOfMassZ0 /= size;
        centerOfMassX1 /= size;
        centerOfMassY1 /= size;
        centerOfMassZ1 /= size;

        centerOfMass0.setXYZ(centerOfMassX0, centerOfMassY0, centerOfMassZ0);
        centerOfMass1.setXYZ(centerOfMassX1, centerOfMassY1, centerOfMassZ1);

        Overlap overlap = new Overlap(coordinates0, coordinates1, size, name0, name1);
        rms = overlap.rms();
        rotateMatrix = overlap.rotationMatrix();

        System.out.println("RMS is " + rms);
        System.out.println("rotation matrix:");
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                System.out.print(rotateMatrix[i][j] + "\t  ");
            }
            System.out.println();
        }

    }


}
