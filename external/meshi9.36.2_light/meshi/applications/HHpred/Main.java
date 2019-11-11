package meshi.applications.HHpred;

import java.io.File;
import java.io.IOException;

import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.AlignmentCell;

public class Main {

	/**
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		String textfile = args[0]; // hhpred file
		String pdb = args[1]; // pdb file
		File file = new File(textfile);
		File pdbfile = new File(pdb);
		HhpredAlignment hhpredalignment = new HhpredAlignment(file);
		AtomList atoms = new AtomList(pdbfile.getAbsolutePath());
		new MolecularSystem();
		Protein protein = new Protein(atoms,
				ResidueExtendedAtomsCreator.creator);
		for (int i = 0; i < hhpredalignment.size(); i++) {
			AlignmentCell tCell = hhpredalignment.get(i).cell(12);
			System.out.println(tCell);
			for (int j = 0; j < protein.atoms().size(); j++) {
				Atom atom = protein.atoms().atomAt(j);
				if (atom.residueNumber() < tCell.number)
					continue;
				else if (atom.residueNumber() > tCell.number)
					break;
				else if (!atom.nowhere()) {
					System.out.println("      " + atom.name + " "
							+ atom.residueName() + " " + atom.residueNumber()
							+ " " + atom.x() + " " + atom.y() + " " + atom.z());
				}
			}

		}

		/*
		 * Writer writer = null; try { writer = new BufferedWriter(new
		 * OutputStreamWriter( new
		 * FileOutputStream("C:\\Users\\Tanya\\Desktop\\RESULTS.txt"),
		 * "utf-8")); writer.write(hal.toString()); } catch (IOException ex) {
		 * // report } finally { try {writer.close();} catch (Exception ex) {} }
		 * bw.write(content.replaceAll("\n",
		 * System.getProperty("line.separator"))); bw.close();
		 */

	}

}
