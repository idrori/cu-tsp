/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.PDB;

import meshi.util.filters.*;

public class PdbLineSEQRES extends PdbLineFilter {
    public boolean acceptPdbLine(PdbLine line) {
        return line.isSEQRES();
    }
}
