package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.ResidueType;
import meshi.util.Utils;
import meshi.util.file.MeshiWriter;

import java.io.File;
import java.io.IOException;
//01234567890123456789
//ATOM    277  CA BSER A  41      18.553 -12.002  40.829  0.50 14.42           C
/**
 * Created by chen on 29/08/2017.
 */
public class FindHolsAndTrpZips {
    public static void main(String[] args) throws IOException{
        MeshiWriter writer = new MeshiWriter("TrpZip.txt");
        MeshiWriter log = new MeshiWriter("log.txt");

        File workingDir = new File(args[0]);
        File[] files = workingDir.listFiles();
        int count = 0;
        for (File file : files) {
            if (count%10 == 0) writer.flush();
            count++;
            if (!file.getName().endsWith(".pdb")) continue;
            System.out.println("Now working on "+file.getName() + " " + count);
            AtomList atoms = new AtomList(file.getName());
            Protein protein = new Protein(atoms,ResidueExtendedAtoms.creator);
            if (atoms.get(0).residue().number() <= 0) {
                log.println("Negative or zero residue number "+file.getName());
            }
            int length;
            ResidueList residues;
            int oneResidueGaps = 0;
            if (protein.chains().size() != 1)  {
                String name = file.getName();
                //file.renameTo(new File(name+".old"));
                MeshiWriter newFile = new MeshiWriter(name+".new");
                char chainChar = name.charAt(4);
                Chain chain = null;
                for (Chain iChain : protein.chains())
                    if (iChain.name().charAt(0) == chainChar) chain = iChain;
                if (chain == null)
                    throw new RuntimeException("This is weird.");
                length = chain.numberOfNonDummyResidues();
                residues = chain.getNonDummyResidues();
                for (Residue residue : chain) {
                    if (!residue.dummy()) {
                        for (Atom atom : residue.getAtoms()) {
                            if (!atom.nowhere())
                                newFile.println(atom);
                        }
                    }
                }
                newFile.close();
               // if (1 == 1) throw new RuntimeException("xxxxxxxxxxxxxxx");
            }
            else {
                length = protein.chain().numberOfNonDummyResidues();
                residues = protein.chain().getNonDummyResidues();
            }
            log.println(file.getName()  + " "+length);

            if (residues.get(residues.size()-1).number()> 995)
                log.println(file.getName()+" top residue number "+residues.get(residues.size()-1).number());
            Residue prevResidue = null;
            for (int iResidue = 0; iResidue < residues.size(); iResidue++ ) {
                Residue residue = residues.get(iResidue);
//                if (residue.type == ResidueType.TRP) {
//                    Atom iCg = residue.cg();
//                    if (!iCg.nowhere()) {
//                        for (int jResidue = iResidue + 1; jResidue < residues.size(); jResidue++) {
//                            Residue secondResidue = residues.get(jResidue);
//                            if (!secondResidue.type().equals(ResidueType.TRP)) continue;
//                            Atom jCg = secondResidue.cg();
//                            if ((!jCg.nowhere()) && iCg.distanceFrom(jCg) < 6.0)
//                                writer.println(protein.name() + " " + residue + " " + secondResidue + " " + iCg.distanceFrom(jCg));
//                        }
//                    }
//                }
                if (prevResidue != null) {
                    if (prevResidue.number() + 1 < residue.number()) {
                        double distance = prevResidue.ca().distanceFrom(residue.ca());
                        String flag = "";
                        if (prevResidue.number() + 2 == residue.number())  {
                            flag = "ONE_RESIDUE_GAP";
                            oneResidueGaps++;
                        }
                        if (prevResidue.number() - residue.number() >= 0)  flag = "WEIRD";
                        log.println("Gap "+protein.name()+" "+prevResidue+" "+residue+" "+distance + " "+ (residue.number() - prevResidue.number()) + " " + flag);
                        if (distance <4.1) {
                            log.println("Gap_with_wierd_distance "+protein.name()+" "+prevResidue+" "+residue+" "+distance);
                        }
                    }
                }
                prevResidue = residue;
            }
            if (oneResidueGaps > 1) log.println(file.getName()+" "+"more then one one-residue-gap");
        }
        writer.close();
        log.close();
    }
}
