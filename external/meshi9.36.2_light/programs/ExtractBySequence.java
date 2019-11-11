package programs;

import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.sequences.*;
import meshi.util.MeshiAttribute;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;
import java.util.Dictionary;
import java.util.Hashtable;
import java.util.Map;

public class ExtractBySequence {
    public static MeshiWriter log;
    public static void main(String[] args) throws IOException{
        Hashtable<String, String> dictionary = new Hashtable<String, String>();
        MeshiLineReader reader = new MeshiLineReader("cb513_pdbs_dict.dat");
        String line;
        while ((line = reader.readLine())!= null) {
            String[] words = line.split(" ");
            dictionary.put(words[1], words[0]);
        }
        reader.close();
        log = new MeshiWriter("extract.log");
        reader = new MeshiLineReader(args[0]);
        while ((line = reader.readLine())!= null) {
            String ID = line.substring(1);
            String name = dictionary.get(ID);
            System.out.println(name+" "+line.substring(1));
            String sequence = reader.readLine();
            reader.readLine();
            extract(name, ID, sequence);
        }
        reader.close();
    }
    public static void extract(String name, String ID, String sequence) throws IOException{
        System.out.println(name+" "+sequence);
        String fileName = name+".pdb";
        String inFile = "pdb\\"+fileName;
        String outFile = "pdbOut\\"+name+"."+ID+".pdb";
        Protein protein = new Protein(new AtomList(inFile), ResidueExtendedAtomsCreator.creator);
        MeshiWriter writer = new MeshiWriter(outFile);

        MeshiSequence meshiSequence = new MeshiSequence(sequence,"");
        new MolecularSystem();
        ResidueSequence proteinSequence = protein.chain().sequence();
        SequenceAlignment sequenceAlignment = SequenceAlignment.identityAlignment(meshiSequence, proteinSequence);
        log.println(name+"\n"+    sequenceAlignment+"\n");
        for (SequenceAlignmentColumn column : sequenceAlignment) {
            AlignmentCell cell = column.cell1();
            Residue residue = (Residue) cell.getAttribute(MeshiAttribute.RESIDUE_ATTRIBUTE);
            if (residue != null)
                residue.getAtoms().print(writer);
        }
        writer.close();
    }
}
