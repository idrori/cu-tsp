package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.parameters.SecondaryStructure;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;

import java.io.IOException;

/**
 */
public class HelixBetaStataistics {
    public static void main(String[] args) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(args[0]);

        String fileName;


        Protein protein;
        int i = 0;
        double minHH = Double.MAX_VALUE, minHO = Double.MAX_VALUE, minOH = Double.MAX_VALUE, minOO = Double.MAX_VALUE;
        double dx, dy, dz, d;
        Atom atomA1, atomA2;
        Atom atomB1, atomB2;
        double dA, dB;
        Residue residueA1, residueA2, residueB1, residueB2;
        int iFile = 0;
        while ((fileName = reader.readLine()) != null) {
            atomA1 = atomA2 = atomB1 = atomB2 = null;
            residueA1 = residueA2 = residueB1 = residueB2 = null;
            dA = dB = Double.MAX_VALUE;
            new MolecularSystem();
            protein = new Protein(new AtomList(fileName), new ResidueExtendedAtomsCreator());
            String dsspFileName = "dssp/" + protein.name() + ".dssp";
            Utils.AssignDSSP(protein, dsspFileName);
            System.out.print(".");
            if (iFile % 200 == 0) System.out.println();
            iFile++;
            for (Atom atom1 : protein.atoms()) {
                if (atom1.nowhere()) continue;
                if ((!atom1.type().backboneH()) && (!atom1.type().backboneO())) continue;
                if (atom1.residue().getSecondaryStructure() != SecondaryStructure.HELIX) continue;
                for (Atom atom2 : protein.atoms()) {
                    if ((!atom2.type().backboneH()) && (!atom2.type().backboneO())) continue;
                    if (atom2.residue().getSecondaryStructure() == SecondaryStructure.HELIX) continue;
                    if ((atom2.residue().getSecondaryStructure() == SecondaryStructure.COIL) && (Math.abs(atom1.residue().number() - atom2.residue().number()) <= 15))
                        continue;
                    if ((atom2.residue().getSecondaryStructure() == SecondaryStructure.SHEET) && (Math.abs(atom1.residue().number() - atom2.residue().number()) <= 4))
                        continue;
                    if (atom2.nowhere()) continue;
                    dx = atom1.x() - atom2.x();
                    dy = atom1.y() - atom2.y();
                    dz = atom1.z() - atom2.z();
                    d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    if (d < 4) {
                        if (atomA1 == null) {
                            atomA1 = atom1;
                            atomA2 = atom2;
                            dA = d;
                            residueA1 = atomA1.residue();
                            residueA2 = atomA2.residue();
                        } else {
                            atomB1 = atom1;
                            atomB2 = atom2;
                            residueB1 = atomB1.residue();
                            residueB2 = atomB2.residue();
                            if (residueA1 != residueB1 &&
                                    residueA2 != residueB2 &&
                                    onTheSameSegment(protein.chain(), residueA1, residueB1) &&
                                    onTheSameSegment(protein.chain(), residueA2, residueB2)) {
                                dB = d;
                                System.out.println("\n_______________________________________\n" +
                                        protein.name() + "\n" +
                                        atomA1 + " " + atomA1.residue().getSecondaryStructure() + "\n" +
                                        atomA2 + " " + atomA2.residue().getSecondaryStructure() + "\n" +
                                        atomB1 + " " + atomB1.residue().getSecondaryStructure() + "\n" +
                                        atomB2 + " " + atomB2.residue().getSecondaryStructure() + "\n" +
                                        dA + "   \t\t    " + dB +
                                        "\n-------------------------------------------------------");
                            }
                        }
                    }
                }
            }
        }

    }

    public static boolean onTheSameSegment(Chain chain, Residue residue1, Residue residue2) {
        if (residue1.number() > residue2.number()) return onTheSameSegment(chain, residue2, residue1);
        SecondaryStructure ss = residue1.getSecondaryStructure();
        int size = chain.size();
        int n1 = residue1.number();
        int n2 = residue2.number();
        for (int i = n1 + 1; (i < size) && (i <= n2); i++) {
            if (chain.get(i).getSecondaryStructure() != ss) return false;
        }
        return true;
    }

}
