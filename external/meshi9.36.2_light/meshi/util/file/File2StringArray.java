/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.file;

import java.io.*;

public class File2StringArray {
    public static String[] f2a(String fileName) {
        String[] result;
        int c = 0;
        try {
            // first pass isOn the file - to find the number of lines
            BufferedReader br = new BufferedReader(new FileReader(fileName));
            String line = br.readLine();
            while (line != null) {
                c++;
                line = br.readLine();
            }
            br.close();
            result = new String[c];
            c = 0;
            // second pass isOn the file - reading the new lines
            br = new BufferedReader(new FileReader(fileName));
            line = br.readLine();
            while (line != null) {
                result[c] = line;
                c++;
                line = br.readLine();
            }
            br.close();
        }
        catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
        return result;
    }
}