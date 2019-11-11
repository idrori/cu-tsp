/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.file;

import java.io.*;
import java.util.*;

import meshi.util.*;

public class MeshiWriter extends PrintWriter {
    private String path = "unknown";

    public MeshiWriter(String fileName) throws IOException {
        super(new BufferedWriter(new FileWriter(fileName)));
    }
    public MeshiWriter(File file) throws IOException {
            super(new BufferedWriter(new FileWriter(file)));
        }

    public MeshiWriter(CommandList commands, Key patheKey, Key nameKey) throws Exception {
        super(getFile(commands, patheKey, nameKey));
        path = (getFile(commands, patheKey, nameKey)).getAbsolutePath();
    }

    public MeshiWriter(CommandList commands, Key patheKey, Key nameKey, String postFix) throws Exception {
        super(getFile(commands, patheKey, nameKey, postFix));
        path = (getFile(commands, patheKey, nameKey, postFix)).getAbsolutePath();
    }


    private static File getFile(CommandList commands, Key patheKey, Key nameKey) {
        return getFile(commands, patheKey, nameKey, "");
    }

    private static File getFile(CommandList commands, Key patheKey, Key nameKey, String postFix) {
        return new File(commands.firstWord(patheKey).secondWord(),
                commands.firstWord(nameKey).secondWord() + postFix);
    }

    public String toString() {
        return path;
    }
}

