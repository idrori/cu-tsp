/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.twoTorsions;

import meshi.energy.Parameters;
import meshi.energy.simpleEnergyTerms.ParametersList;
import meshi.geometry.Torsion;
import meshi.geometry.TorsionPair;
import meshi.molecularElements.Residue;

import java.util.Iterator;

public class TwoTorsionsParametersList extends ParametersList {

    public TwoTorsionsParametersList(String[] parametersFileName) {
        super(parametersFileName, false);  // non-sortable
    }

    public Parameters createParameters(String line) {
        return new TwoTorsionsParameters(line);
    }

    public Parameters parameters(Object obj) {
        TorsionPair torsionPair = (TorsionPair) obj;
        Torsion torsion1 = torsionPair.torsion1();
        Torsion torsion2 = torsionPair.torsion2();
        Residue residue = torsion1.atom2.residue();
        int resnum1 = torsion1.getTorsionResNum();
        int resnum2 = torsion2.getTorsionResNum();
        int code1 = torsion1.getTorsionCode();
        int code2 = torsion2.getTorsionCode();
        String name1 = torsion1.getTorsionName();
        String name2 = torsion2.getTorsionName();
        TwoTorsionsParameters twoTorsionsParameters;
        if (resnum1 != resnum2) return null;
        if (code1 < 0) return null;
        if (code2 < 0) return null;
        Iterator iter = iterator();
        while (iter.hasNext()) {
            twoTorsionsParameters = (TwoTorsionsParameters) iter.next();
            if ((name1.equals(twoTorsionsParameters.torsion1Name)) &&
                    (name2.equals(twoTorsionsParameters.torsion2Name)) &&
                    (twoTorsionsParameters.mapping[residue.type.ordinal()] > -1) &&
                    (!(name2.equals("PHI") || name2.equals("PSI") || name2.equals("CHI1")) ||
                            residue.getSecondaryStructure() == twoTorsionsParameters.secondaryStructure || // Secondary structure exact match
                            residue.getSecondaryStructure().equalsIgnorOr(twoTorsionsParameters.secondaryStructure))) // XXX_OR_COIL match to COIL
                return twoTorsionsParameters;
        }
        return null;
    }
}
