/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.file;

import meshi.util.*;
import meshi.util.string.*;
import meshi.util.filters.*;

import java.io.*;
import java.util.*;
import java.util.zip.*;
// The major tool for reading asci and gziped files in meshi
// #1 variables
// #2 constructors

// #3 methods
//  #3a constructor supplements
//  #3b path 
public class MeshiLineReader extends LineNumberReader {
    // #1 variables
    private String path; // the path to the file
    private String name;
    private File file = null;

    public File file() {
        return file;
    }

    // #2 constructors
    // The different constructors are used to funnel the different 
    // java classes for input reading into a single class.
    public MeshiLineReader(InputStreamReader ISR, String path) {
        super(ISR);
        setPath(path);
        file = new File(path);
    }

    public MeshiLineReader(InputStreamReader ISR) {
        this(ISR, "Unknown path");
    }

    public MeshiLineReader(FileReader FR, String path) {
        super(FR);
        setPath(path);
        file = new File(path);
    }

    public MeshiLineReader(FileReader FR) {
        this(FR, "Unknown path");
    }

    public MeshiLineReader(String path) {
        this(new File(path));
        file = new File(path);
    }

    // Note that this constructor allows us to read gzipped files
    // "without noticing it".
    public MeshiLineReader(File file) {
        this((isItGzipped(file) ?
                new_InputStreamReader(file) :
                new_FileReader(file)),
                file.getAbsolutePath());
        this.file = file;
    }

    // #3a constructor supplements
    private static boolean isItGzipped(File file) {
        return file.getAbsolutePath().endsWith(".gz");
    }

    private static InputStreamReader new_InputStreamReader(File file) {
        try {
            return new InputStreamReader(new GZIPInputStream(new FileInputStream(file)));
        }
        catch (Exception e) {
            throw new MeshiException("MeshiLineReader Error while trying to open gziped file:" +
                    "File name: " + file.getAbsolutePath() + "\n" +
                    e.getMessage() + "\n");
        }
    }

    private static FileReader new_FileReader(File file) {
        try {
            return new FileReader(file);
        }
        catch (Exception e) {
            throw new MeshiException("MeshiLineReader Error while trying to open file:" +
                    "File name: " + file.getAbsolutePath() + "\n" +
                    e.getMessage() + "\n");
        }
    }

    // #3b path
    private void setPath(String p) {
        path = p;
        StringList separators = StringParser.bySeparator("/ \\ ", " ");
        StringList pathAsAList = StringParser.bySeparators(path, separators);
        name = pathAsAList.get(pathAsAList.size() - 1);
    }

    public String path() {
        return path;
    }

    public String fileName() {
        StringList temp = StringParser.bySeparator(path, "/");
        return temp.get(temp.size() - 1);
    }

    public String name() {
        return name;
    }

    public String readLine(Filter filter) {
        String temp;
        try {
            while ((temp = readLine()) != null)
                if (filter.accept(temp)) return temp;
        }
        catch (Exception e) {
            throw new MeshiException("MeshiLineReader Error while trying to read file:" +
                    "File name: " + path() + "\n" +
                    e.getMessage() + "\n");
        }
        return null;

    }

    public String readLine(String commentString) {
        String temp;
        try {
            while ((temp = readLine()) != null) {
                if (!temp.startsWith(commentString)) {
                    if (temp.indexOf(commentString) == -1) return temp;
                    else return temp.substring(0, temp.indexOf(commentString));
                }
            }
        }
        catch (Exception e) {
            throw new MeshiException("MeshiLineReader Error while trying to read file:" +
                    "File name: " + path() + "\n" +
                    e.getMessage() + "\n");
        }
        return null;

    }

    public String toString() {
        return "MeshiLineReader";
    }
}
