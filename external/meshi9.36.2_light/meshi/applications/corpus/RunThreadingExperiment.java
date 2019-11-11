/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.corpus;

import meshi.energy.EnergyCreator;
import meshi.util.*;

public class RunThreadingExperiment extends MeshiProgram {

    private static String corpusFileName = null;

    private static int Nexp = -1;

    private static int Ninst = -1;

    private static int fragL = -1;

    private static int RMSind1 = -1;

    private static int RMSind2 = -1;


    public static void main(String[] args) {
        init(args);
        Corpus corpus = new Corpus(corpusFileName);
        corpus.findPhiPsiRelations(100000);
        if (1 == 1)
            System.exit(0);
        for (int c = 0; c < Nexp; c++)
            corpus.GeneralThreadingExperiment(fragL, Ninst, c + "", RMSind1, RMSind2);
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
                "Usage java -Xmx800m RunThreadingExperiment <random seed> <corpus filename> <Num of Experiments> <Num of instances in experiment> <fragL> <RMSind1> <RMSind2>\n" +
                "                    ******************\n");

        if (getFlag("-debug", args)) tableSet("debug", new Boolean(true));

        String seedString = getOrderedArgument(args);
        if (seedString == null) throw new RuntimeException(errorMessage);
        int seed = (new Integer(seedString)).intValue();
        System.out.println("# seed is " + seed);
        initRandom(seed);

        corpusFileName = getOrderedArgument(args);
        if (corpusFileName == null) throw new RuntimeException(errorMessage);
        System.out.println("# Corpus file name is " + corpusFileName);

        line = getOrderedArgument(args);
        if (line == null) throw new RuntimeException(errorMessage);
        Nexp = (new Integer(line)).intValue();
        System.out.println("# Nexp " + Nexp);

        line = getOrderedArgument(args);
        if (line == null) throw new RuntimeException(errorMessage);
        Ninst = (new Integer(line)).intValue();
        System.out.println("# Ninst " + Ninst);

        line = getOrderedArgument(args);
        if (line == null) throw new RuntimeException(errorMessage);
        fragL = (new Integer(line)).intValue();
        System.out.println("# fragL " + fragL);

        line = getOrderedArgument(args);
        if (line == null) throw new RuntimeException(errorMessage);
        RMSind1 = (new Integer(line)).intValue();
        System.out.println("# RMSind1 " + RMSind1);

        line = getOrderedArgument(args);
        if (line == null) throw new RuntimeException(errorMessage);
        RMSind2 = (new Integer(line)).intValue();
        System.out.println("# RMSind2 " + RMSind2);
    }
}
