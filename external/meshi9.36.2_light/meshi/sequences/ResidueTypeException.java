/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences;

public class ResidueTypeException extends RuntimeException {
    public ResidueTypeException(Character weirdChar, MeshiSequence sequence) {
        super("Weird character " + weirdChar + " in \n" + sequence);
    }

    public ResidueTypeException(Character weirdChar, String sequence) {
        super("Weird character " + weirdChar + " in \n" + sequence);
    }

    public ResidueTypeException(Character weirdChar) {
        super("Weird character " + weirdChar);
    }
}
