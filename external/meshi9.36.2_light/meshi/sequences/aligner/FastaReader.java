/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 28/07/2005
 * Time: 13:23:52
 * To change this template use File | Settings | File Templates.
 */

import java.io.*;
import java.util.*;

/**
 * Text Files: read, write:<br><br>
 * this class covers some importent topics of cs such as text-files: <br>
 * Note: this code contains a wrong exception handling!!
 */
public class FastaReader {
    String str = new String();

    public FastaReader() {
    }

    public FastaReader(String file) {

        try {
            FileReader fr = new FileReader(file);
            BufferedReader is = new BufferedReader(fr);
            String s1 = is.readLine();
            while (s1 != null) {
                str = str.concat("\n" + s1);
                s1 = is.readLine();
            }
            is.close();


        }

        catch (IOException e) {
            System.out.print("Error: while reading" + e);
            System.exit(2);
        }
        // System.out.println(fstr.getSource());
        //       System.out.println(fstr.getImage());
    }

    public String getSource() {
        String src = new String();
        StringTokenizer tokenizer = new StringTokenizer(str, ">", false);
        if (tokenizer.hasMoreElements())
            src = tokenizer.nextToken();  //gets  first > which after that will be the sequence were looking for
        else System.out.println("wrong input1");//TODO change to exceptions
        if (tokenizer.hasMoreElements()) src = tokenizer.nextToken();
        else System.out.println("wrong input1");
        tokenizer = new StringTokenizer(src, "\n", false);
        if (tokenizer.hasMoreElements()) src = tokenizer.nextToken();
        else System.out.println("wrong input2");//TODO change to exceptions
        if (tokenizer.hasMoreElements()) src = tokenizer.nextToken();
        else System.out.println("wrong input2");//TODO change to exceptions
        while (tokenizer.hasMoreElements())
            src = src + tokenizer.nextToken();
        return src;
    }

    public String getImage() {
        String img = new String();
        StringTokenizer tokenizer = new StringTokenizer(str, ">", false);
        if (tokenizer.hasMoreElements())
            img = tokenizer.nextToken();  //gets  first > which after that will be the sequence were looking for
        else System.out.println("wrong input1");//TODO change to exceptions
        if (tokenizer.hasMoreElements()) img = tokenizer.nextToken();
        else System.out.println("wrong input1");
        if (tokenizer.hasMoreElements()) img = tokenizer.nextToken();
        else System.out.println("wrong input1");
        tokenizer = new StringTokenizer(img, "\n", false);
        if (tokenizer.hasMoreElements()) img = tokenizer.nextToken();
        else System.out.println("wrong input2");//TODO change to exceptions
        if (tokenizer.hasMoreElements()) img = tokenizer.nextToken();
        else System.out.println("wrong input2");//TODO change to exceptions
        while (tokenizer.hasMoreElements())
            img = img + tokenizer.nextToken();
        return img;
    }
    //  public static void  main (String [] a)    {
    // FastaReader fstr=new FastaReader("Q9SZR1.fas") ;
    //    System.out.println(fstr.getSource());
    //    System.out.println(fstr.getImage());

//}
}