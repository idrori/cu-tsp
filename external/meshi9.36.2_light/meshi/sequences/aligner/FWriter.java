/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 07/08/2005
 * Time: 20:46:46
 * To change this template use File | Settings | File Templates.
 */

import java.io.IOException;
import java.io.PrintWriter;
import java.io.FileWriter;


public class FWriter {
    PrintWriter os;

    public FWriter(String out) {
        try {
            FileWriter fi = new FileWriter(out);
            os = new PrintWriter(fi);
        }
        catch (IOException e) {
            System.out.print("Error while writing: " + e);
            System.exit(1);
        }
    }

    public void write(String word) {
        os.print(word);
    }

    public void write(char word) {
        os.print(word);
    }

    public void writeln() {
        os.println();
    }

    public void writeln(String line) {
        os.println(line);
    }

    public void close() {
        os.close();
    }
}