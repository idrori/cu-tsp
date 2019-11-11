/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.*;
import meshi.util.*;
import meshi.PDB.*;

import java.util.*;

import meshi.energy.simpleEnergyTerms.bond.*;
import meshi.energy.simpleEnergyTerms.angle.*;

/**
 * An extended getAtoms model of a protein.
 * This model includes all heavy getAtoms and most hydrogen getAtoms that take part in hydrogen bonds.
 */
public class ExtendedAtomsProtein extends Protein {
    /**
     * The parameter is expected to be a PDB formatted file.
     */
    public ExtendedAtomsProtein(String fileName, CommandList commands) {
        this(fileName, commands, new ArrayList<String>());
    }

    public ExtendedAtomsProtein(String fileName, CommandList commands, List<String> chainNames) {
        super(fileName, new PdbLineATOM(chainNames), new ResidueExtendedAtomsCreator(), null);

        BondParametersList bondParametersList = Utils.getBondParameters(commands);
        AngleParametersList angleParametersList = Utils.getAngleParameters(commands);
        for (Iterator iter = residues().iterator(); iter.hasNext();) {
            Residue residue = (Residue) iter.next();
            if (!residue.dummy())
                ((ResidueExtendedAtoms) residue).addHydrogens(bondParametersList, angleParametersList);
        }
    }

}
