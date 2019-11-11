/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 */
public class MeshiInfoXMLwriter extends MeshiWriter {
    public static final String BLANK10 = "          ";
    public static final String BLANK20 = "                    ";
    public static final String BLANK30 = "                              ";
    public MeshiInfoXMLwriter(String fileName) throws IOException {
        super(fileName);
        printf("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
    }
}
