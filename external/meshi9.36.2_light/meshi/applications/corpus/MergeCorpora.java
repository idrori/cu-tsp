/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import meshi.parameters.*;
import meshi.util.*;
import meshi.util.file.File2StringArray;

public class MergeCorpora extends MeshiProgram { /**
 * The implemented
 * interfaces defines the
 * names of atom and residue
 * types.
 **/

    /**
     * A string with the name of the pdb file to minimize.
     */
    private static String fileName = null;


    public static void main(String[] args) {
        init(args);
        String[] files = File2StringArray.f2a(fileName);
        Corpus corpus = new Corpus(files[0]);
        for (int c = 1; c < files.length; c++) {
            System.out.println("Merging: " + files[c]);
            Corpus tmpCorpus = new Corpus(files[c]);
            corpus.merge(tmpCorpus);
        }
        corpus.writeToDisk("finalCorpus.txt");
        //corpus.threadingExperiment(10);
        //Corpus corpus = new Corpus(fileName);
        //corpus.threadingExperiment_withPP(10,1.0,0.0);
    }

    /**
     * ================================= init =========================================
     * <p/>
     * A static function for parsing of the command line arguments and assigning the
     * variables commandsFileName, modelFileName and randomNumberSeed with the right inputs. Note that this
     * static method is using parsing functions such as getOrderedArguments that are defined in MeshiProgram
     * that MinimizeProtein inherits.
     */

    protected static void init(String[] args) {


        String line;
        String errorMessage = ("\n                  ******************\n" +
                "Usage java -Xmx300m MergeCorpora <corpora file> \n" +
                "                    ******************\n");

        if (getFlag("-debug", args)) tableSet("debug", new Boolean(true));

        fileName = getOrderedArgument(args);
        if (fileName == null) throw new RuntimeException(errorMessage);
        System.out.println("# initial model file name is " + fileName);

        initRandom(333);
    }
}
