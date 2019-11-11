package programs;

import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by chen on 29/11/2015.
 */
public class extractChainsFromPDBfiles {
    public static void main(String[] args) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(args[0]);
        String line;
        reader.readLine();
        ArrayList<NameLengthPair> list = new ArrayList();
        while ((line = reader.readLine()) != null) {
            String[] words = line.split(" +");
            NameLengthPair nameLengthPair = new NameLengthPair(words[0].substring(0, 4), words[0].substring(4, 5), Integer.valueOf(words[1]).intValue());
            list.add(nameLengthPair);
        }
        reader.close();
        Object[] array = list.toArray();
        Arrays.sort(array);
        for (Object o : array) {
            System.out.println(o);
            NameLengthPair nameLengthPair = (NameLengthPair) o;
            String outFileName = "..\\chains\\"+nameLengthPair.name.substring(0,4)+nameLengthPair.chain+".pdb";
            File file = new File(outFileName);
            if (! file.exists()) {
                AtomList atomList = new AtomList(nameLengthPair.name);
                Protein protein = new Protein(atomList, ResidueExtendedAtomsCreator.creator);
                Chain chain = protein.getChain(nameLengthPair.chain);
                MeshiWriter writer = new MeshiWriter("..\\chains\\" + nameLengthPair.name.substring(0, 4) + nameLengthPair.chain + ".pdb");
                chain.atoms().print(writer);
                writer.close();
            }
        }
    }

//            String proteinName = words[0].substring(0, 4) + ".pdb";
//            File file = new File(proteinName);
//            if (!file.exists()) {
//                URL url = new URL(" \thttp://www.rcsb.org/pdb/files/" + proteinName);
//                System.out.println(url.toString());
//                InputStream inputStream = url.openStream();
//                InputStreamReader inputStreamReader = new InputStreamReader(inputStream);
//                BufferedReader bufferedReader = new BufferedReader(inputStreamReader);
//                MeshiWriter writer = new MeshiWriter(proteinName);
//                while ((line = bufferedReader.readLine()) != null)
//                    writer.println(line);
//                writer.close();
//                bufferedReader.close();
//                inputStreamReader.close();
//                inputStream.close();
//            }
//        }
//    }

    private static class NameLengthPair implements Comparable <NameLengthPair> {
        String name;
        String chain;
        int length;
        public NameLengthPair(String name, String chain, int length ){
            this.chain = chain;
            this.length = length;
            this.name = name;
        }

        public int compareTo(NameLengthPair nameLengthPair) {
            return length - nameLengthPair.length;
        }

        public String toString() {
            return name+" "+chain+" "+length;
        }
    }
}
