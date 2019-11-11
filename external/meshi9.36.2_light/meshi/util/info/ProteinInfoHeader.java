/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.molecularElements.Protein;

/**
 *
 */
public class ProteinInfoHeader extends StringInfoElement {
    private String fileName;
    private String name;
    private Protein protein;

    public ProteinInfoHeader(String name, String fileName, String comment, Protein protein) {
        super(InfoType.ProteinInfo, comment, name);
        this.fileName = fileName;
        this.name = name;
        this.protein = protein;
    }

    public String toXML() {
        String out = "<" + type +
                                  "  value=\"" + name + "\" "+
                                  "fileName=\"" + fileName + "\" "+
                                 // "length=\"" + protein.chain().numberOfNonDummyResidues() + "\" "+
                                  ">";
        return out;
    }

}
