package meshi.applications.HHpred;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;

public class tiuta {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String pdb = args[0]; // pdb file

		String content = "";
		File pdbfile = new File(pdb);
		AtomList atoms = new AtomList(pdbfile.getAbsolutePath());
		new MolecularSystem();
		Protein protein = new Protein(atoms,
				ResidueExtendedAtomsCreator.creator);
		// System.out.println(protein.getAtoms().toString());
		for (Residue res : protein.residues()) {
			content = content + res.name + "  " + res.number() + "\n";
			System.out.println(res.name + "  " + res.number());
			for (Atom atom : res.getAtoms()) {
				content = content + atom.toString() + "\n";
				System.out.println(atom.toString());
			}
		}
		Writer writer = null;
		try {
			writer = new BufferedWriter(new OutputStreamWriter(
					new FileOutputStream(
							"C:\\Users\\Tanya\\Desktop\\protein.txt"), "utf-8"));
			writer.write(content.replaceAll("\n",
					System.getProperty("line.separator")));
		} catch (IOException ex) {
			// report
		} finally {
			try {
				writer.close();
			} catch (Exception ex) {
			}
		}

	}

}
