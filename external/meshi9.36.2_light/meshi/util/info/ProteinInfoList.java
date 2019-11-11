/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.ArrayList;

/**

 */
public class ProteinInfoList extends ArrayList<ProteinInfo> {
    public String toXml(int indentation, String name) throws IOException {
        String out = indentation+"<ProteinInfoList name=\""+name+"\">\n";
        for (ProteinInfo proteinInfo : this) {
            out += proteinInfo.toXml(indentation+1);
        }
        out += "</ProteinInfoList>";
        return out;
    }

    public void writeXml(MeshiWriter writer, String name) throws IOException {
        writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>");
        writer.println(toXml(3, name));
    }
}
