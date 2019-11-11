/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.util.formats.Format;
import meshi.util.string.StringList;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import java.util.StringTokenizer;

/**
 * The (suggested) superclass of programs using the meshi package.
 * It includes a bunch of usefull static methods and global variables.
 */
public class MeshiProgram {
    //----------------------------- class variables ---------------------------------------------
    /**
     * After program initiation, stores all the information needed for program execution.
     * It is arranged as a list of (key,value) pairs
     */
    private static GlobalElementList globalTable = new GlobalElementList();
    /**
     * Program name.
     */
    protected static String name = "MeshiProgram";
    /**
     * The command line with which the program was called.
     */
    protected static String commandLine = null;
    /**
     * The single random number generator used by the program.
     */
    private static Random randomNumberGenerator = null;
    /**
     * The seed of the  random number generator.
     */
    private static int seed;

    public static int seed() {
        return seed;
    }

    // ------------------------------- command Line parsing --------------------------------------

    /**
     * Get a boolean flag
     * Searcheds for a the keyword <b>key</b> in <b>args</b>.
     * If the keyword is found returns true and removed from <b>args</b>.
     * If the keyword is not found in <b> args</b> returns false.
     */
    protected static boolean getFlag(String key, String[] args) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals(key)) {
                args[i] = "";
                if (getFlag(key, args)) {
                    throw new MeshiException("Flag " + key + " appear more then onece in command line:\n" +
                            commandLine + "\n");
                }
                return true;
            }
        }
        return false;
    }

    /**
     * Searches for a the keyword <b>key</b> in <b>args</b>.
     * If the keyword is found, it is removed from <b>args</b>. The following word (the associated value)
     * is also removed and  returned.
     * If the keyword is not found in <b> args</b>null is returned.
     */
    protected static String getFlagedArgument(String key, String[] args) {
        String temp;
        for (int i = 0; i < args.length - 1; i++) {
            if (args[i].equals(key)) {
                args[i] = "";
                temp = args[i + 1];
                args[i + 1] = "";
                if (getFlag(key, args)) {
                    throw new MeshiException("Flag " + key + " appear more then onece in command line:\n" +
                            commandLine + "\n");
                }
                return temp;
            }
        }
        return null;
    }

    public static String getArgument(String key, String[] args) {
        String temp = getArgument(key,args,0);
//        if (temp == null) {
//            System.out.println("Arguments are:");
//            for (String arg : args) {
//                System.out.println(arg);
//            }
//            throw new MeshiException(key + " not found in command line.");
//        }
        return temp;
    }

    protected static String[] getArguments(String[] keys,String[] args){
        String[] out = new String[keys.length];
        try {
            for (int i = 0; i < keys.length; i++) {
                String key = keys[i];
                out[i] = getArgument(key,args);
            }
        }
        catch (MeshiException ex) {
            System.err.print("Required arguments: ");
            for (String key : keys) System.err.print(" "+key);
            System.err.println();
            System.err.print("found\n");
            for (String arg : args) System.err.print(" "+arg);
            System.err.println();
        }

        return out;
    }
    private static String getArgument(String key, String[] args, int start) {
        String temp = null;
        for (int i = start; i < args.length; i++) {
            String fullKey = "-"+key+"=";
            if (args[i].startsWith(fullKey)) {
                temp = args[i].substring(fullKey.length());
                if (i < args.length-1) {
                    if (getArgument(key,args,i+1)!=null)
                            throw new MeshiException("Flag " + key +
                                                     " appear more then once in command line:\n");
                }
                return temp;
            }
        }
        return null;
    }

    /**
     * Removes the first string in <b>args</b> and returns it.
     */
    protected static String getOrderedArgument(String[] args) {
        String temp;
        for (int i = 0; i < args.length; i++) {
            if (!((temp = args[i]).equals(""))) {
                args[i] = "";
                return temp;
            }
        }
        return null;
    }

    //-------------------------------------- setting parameters -----------------------------------

    /**
     * Associate a value to a keyword.
     */
    protected static Object tableSet(String key, Object value) {
        return (globalTable.set(key, value).value());
    }

    //-------------------------------------- getting parameters -----------------------------------

    public static boolean tableIncludes(String key) {
        GlobalElement element;
        element = globalTable.get(key);
        if (element == null) return false;
        else return true;
    }

    /**
     * Gets the value associated with the keyword <b>key</b>.
     */
    protected static Object tableGet(String key) {
        GlobalElement element;
        if (globalTable == null)
            throw new RuntimeException("globalTable not initialized. Its fine but please do not \n" +
                    "try to read it");
        element = globalTable.get(key);
        if (element == null) {
            throw new MeshiException("Global table error:\n" +
                    "No element with key = \"" +
                    key + "\"\n" +
                    globalTable);
        }
        return (element.value());
    }

    /**
     * Gets the Double value associated with the keyword <b>key</b>.
     */
    public static Double getD(String key) {
        try {
            return (Double) tableGet(key);
        }
        catch (Exception e) {
            throw new RuntimeException("Error in getD(" + key + ")\n" + e);
        }
    }

    /**
     * Gets the Integer value associated with the keyword <b>key</b>.
     */
    public static Integer getI(String key) {
        try {
            return (Integer) tableGet(key);
        }
        catch (Exception e) {
            throw new RuntimeException("Error in getI(" + key + ")\n" + e);
        }
    }

    /**
     * Gets the Boolean value associated with the keyword <b>key</b>.
     */
    public static Boolean getB(String key) {
        try {
            return (Boolean) tableGet(key);
        }
        catch (Exception e) {
            throw new RuntimeException("Error in getB(" + key + ")\n" + e);
        }
    }

    /**
     * Gets the String value associated with the keyword <b>key</b>.
     */
    public static String getS(String key) {
        try {
            return (String) tableGet(key);
        }
        catch (Exception e) {
            throw new RuntimeException("Error in getS(" + key + ")\n" + e);
        }
    }

    public static String getS(String key, StringList commands, String commandsFileName) {
        String line;
        line = commands.startsWith(key);
        if (line == null)
            throw new RuntimeException("A problem while reading commands file " + commandsFileName + "\n" +
                    "No line starts with the keyword " + key + ".");
        return line;
    }

    protected static String get2ndString(String line, String commandsFileName) {
        StringTokenizer tokenizer;
        try {
            tokenizer = new StringTokenizer(line);
            tokenizer.nextToken(); // 1st word
            return tokenizer.nextToken();
        }
        catch (Exception e) {
            throw new RuntimeException("A problem while reading commands file " + commandsFileName + "\n" +
                    "cannot parse:\n\t" + line + "\n" + "second string not found\n" + e);
        }
    }


    /**
     * Gets the double value associated with the keyword <b>key</b>.
     */
    public static double getd(String key) {
        return getD(key).doubleValue();
    }

    /**
     * Gets the int value associated with the keyword <b>key</b>.
     */
    public static int geti(String key) {
        return getI(key).intValue();
    }

    /**
     * Gets the boolean value associated with the keyword <b>key</b>.
     */
    public static boolean getb(String key) {
        return getB(key).booleanValue();
    }

    /**
     * Gets the boolean value associated with the keyword <b>verbose</b>.
     */
    public static boolean verbose() {
        if (globalTable == null) return false;
        if (!tableIncludes("verbose")) return false;
        return getb("verbose");
    }

    /**
     * Gets the boolean value associated with the keyword <b>debug</b>.
     */
    public static boolean debug() {
        if (globalTable == null) return false;
        if (!tableIncludes("debug")) return false;
        return getb("debug");
    }

    /**
     * Gets the program's random number generator.
     */
    public static Random randomNumberGenerator() {
        if (randomNumberGenerator == null)
            throw new RuntimeException("The randomNumberGenerator was not initialized use \n" +
                    "initRandom() or initRandom(int seed).\n" +
                    "initRandom() is equivalent to initRandom(0)");
        return randomNumberGenerator;
    }


    //------------------------------------- helper methods ----------------------------------------

    public String toString() {
        if (globalTable == null) return name;
        Iterator iter = globalTable.iterator();
        GlobalElement element;
        String out = name + "\n" + "parameters:\n";
        while ((element = (GlobalElement) iter.next()) != null)
            out += element.key() + "\t" + element.value() + "\n";
        return out;
    }

    /**
     * Initializes the random number generator.
     */
    protected static void initRandom() {
        if (randomNumberGenerator != null)
            throw new RuntimeException("The random number generator may be initialized only once");
        MeshiProgram.seed = seed;
        randomNumberGenerator = new Random();
    }

    public static void initRandom(int seed) {
        if (randomNumberGenerator != null)
            throw new RuntimeException("The random number generator may be initialized only once");
        MeshiProgram.seed = seed;
        randomNumberGenerator = new Random(seed);
    }

    protected static int initRandom(String[] argv) {
        String arg = getFlagedArgument("-seed", argv);
        if (arg == null) seed = 0;
        else seed = (new Integer(arg)).intValue();
        initRandom(seed);
        return seed;
    }

    /**
     * What you always wanted to know.
     */
    public static String about() {
        return "" +
                "The meshi package is developed in Chen Keasars Lab \n" +
                "at the departments of Computer Science and Life Siences,\n" +
                "Ben-Gurion University of the Negev, Be'er-Sheva, Israel\n\n" +
                "Contributors (in alphabetical order):\n" +
                "Arik David\n" +
                "Ohad Givati\n" +
                "Chen Keasar\n" +
                "Danny Klein\n" +
                "Diana Lavie\n" +
                "Lena Margolis\n" +
                "Ofer Meisles\n" +
                "Yael Pinchasov\n" +
                "Elior Siebert\n" +
                "Sharon Zafriri\n" +
                "Ziv Zeira\n" +
                "Tal-Sarit Zobakov";
    }

    /**
     * Assembles the command-line string from an array of command-line-words.
     */
    private static String getCommandLine(String[] args) {
        String commandLine = "";
        for (int i = 0; i < args.length; i++)
            commandLine = commandLine + " " + args[i];
        return commandLine;
    }

    /**
     * Print the keyword - value table.
     */
    public static void printGlobalTable() {
        globalTable.print();
    }
    //----------------------------------------- helper classes ---------------------------------------

    /**
     * An element of the keyword-value table.
     */
    private static class GlobalElement {
        String key;
        Object value;

        public GlobalElement(String key, Object value) {
            this.key = key;
            this.value = value;
        }

        public String key() {
            return key;
        }

        public Object value() {
            return value;
        }

        public void setValue(Object value) {
            this.value = value;
        }
    }

    /**
     * A keyword-value table
     */
    private static class GlobalElementList extends ArrayList {
        public GlobalElementList() {
            if (globalTable != null)
                throw new RuntimeException("Global table may not be created more than onece.");
        }

        private GlobalElement add(String key, Object value) {
            GlobalElement newElement = new GlobalElement(key, value);
            add(newElement);
            return newElement;
        }

        public GlobalElement get(String key) {
            Iterator iter = iterator();
            GlobalElement element;
            while ((element = (GlobalElement) iter.next()) != null)
                if (element.key().equals(key)) return element;
            return null;
        }

        public GlobalElement set(String key, Object value) {
            GlobalElement element = get(key);
            if (element == null) return add(key, value);
            else element.setValue(value);
            return element;
        }

        public String toString() {
            Format sformat = Format.STANDARD;
            Iterator elements = iterator();
            GlobalElement element;
            String out = "Global elements table:\n";
            while ((element = (GlobalElement) elements.next()) != null)
                out += sformat.f(element.key()) + " " + element.value() + "\n";
            return out;
        }

        public void print() {
            System.out.println(toString());
        }
    }
}
