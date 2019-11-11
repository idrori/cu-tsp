/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueCreator;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;
import meshi.parameters.BBatom;
import meshi.parameters.ResidueMode;
import meshi.parameters.ResidueType;

/**
 * <pre>
 * Alanin from Levitt, JMB 168:592 (1983) table 2.
 *           CB  O
 *           |   |
 *      N - CA - C...n
 *      |
 *      H
 */
public class BackboneResidue extends ResidueExtendedAtoms implements ResidueCreator {
    public static final String COMMENT = "General backbone residue.";
    public static final BackboneResidue creator = new BackboneResidue("creator");

    public BackboneResidue(String name) {
        super(name);
    }

    public BackboneResidue(ResidueType type, AtomList atomList,
                           ResidueIdentifier id, ResidueMode mode,
                           MolecularSystem molecularSystem) {
        super(type, atomList, id, mode,molecularSystem);
    }

    public String comment() {
        return COMMENT;
    }

    public Residue create(AtomList atoms, ResidueIdentifier id,
                          ResidueMode mode,MolecularSystem molecularSystem) {
        if ((atoms == null) || (atoms.size() == 0)) return new Residue(id);
        Atom atom = find(atoms, BBatom.CA);
        AtomType atomType = atom.type();
        ResidueType residueType = ResidueType.type(atomType);
        return new BackboneResidue(residueType, atoms, id, mode,molecularSystem);
    }
}
