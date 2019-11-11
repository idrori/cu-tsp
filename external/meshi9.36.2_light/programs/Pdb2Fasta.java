package programs;

import meshi.util.Utils;

import java.io.File;
import java.io.IOException;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 27/06/2010
 * Time: 14:41:23
 * To change this template use File | Settings | File Templates.
 */
public class Pdb2Fasta {
    public static void main(String[] argv) throws IOException{
        File dir = new File(".");
        File[] files = dir.listFiles();
        for (File file : files){
            if (file.getName().endsWith("pdb")) {
                Utils.pdb2Fasta(file);
            }
        }
    }

}
