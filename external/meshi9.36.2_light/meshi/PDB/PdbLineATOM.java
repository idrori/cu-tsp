/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB;

import meshi.util.filters.*;

import java.util.*;

public class PdbLineATOM extends PdbLineFilter {
    public final List<String> chainNames;
    public static final Filter filter = new PdbLineATOM();

    public PdbLineATOM() {
        this(new ArrayList<String>());
    }

    public PdbLineATOM(List<String> chainNames) {
        this.chainNames = chainNames;
    }

    public boolean acceptPdbLine(PdbLine line) {
        if (chainNames.size() == 0) return line.isAnAtom();
        for (String name : chainNames) {
            if (line.isAnAtom() && line.chain().equals(name)) return true;
        }
        return false;
    }
}
