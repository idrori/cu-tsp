/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.extendedAtoms;

import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.geometry.*;
import meshi.util.*;

import java.util.*;

public class ResidueExtendedAtomsCreator implements ResidueCreator {
    public static final ResidueExtendedAtomsCreator creator = new ResidueExtendedAtomsCreator();

    public ResidueExtendedAtomsCreator() {
    }

    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        if (atoms == null)
            throw new RuntimeException("This is weird");
        Atom atom = atoms.atomAt(0);
        String name;
        ResidueType type;
        if (atom == null) throw new RuntimeException("This is weird " + id);
        if (atom.residue() != null) {
            name = atom.residue().name;
            type = ResidueType.type(name);
        } else {
            type = ResidueType.type(atom.type());
        }
        return create(atoms, type, id, mode, molecularSystem);
    }


    public Residue create(ResidueType type, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
         return create(null,type,id,mode,molecularSystem);
    }

    public Residue create(AtomList atoms, ResidueType type, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        if (type == ResidueType.ALA) return new Ala(atoms, id, mode,molecularSystem);
        if (type == ResidueType.CYS) return new Cys(atoms, id, mode,molecularSystem);
        if (type == ResidueType.ASP) return new Asp(atoms, id, mode,molecularSystem);
        if (type == ResidueType.GLU) return new Glu(atoms, id, mode,molecularSystem);
        if (type == ResidueType.PHE) return new Phe(atoms, id, mode,molecularSystem);
        if (type == ResidueType.GLY) return new Gly(atoms, id, mode,molecularSystem);
        if (type == ResidueType.HIS) return new His(atoms, id, mode,molecularSystem);
        if (type == ResidueType.ILE) return new Ile(atoms, id, mode,molecularSystem);
        if (type == ResidueType.LYS) return new Lys(atoms, id, mode,molecularSystem);
        if (type == ResidueType.LEU) return new Leu(atoms, id, mode,molecularSystem);
        if (type == ResidueType.MET) return new Met(atoms, id, mode,molecularSystem);
        if (type == ResidueType.ASN) return new Asn(atoms, id, mode,molecularSystem);
        if (type == ResidueType.PRO) return new Pro(atoms, id, mode,molecularSystem);
        if (type == ResidueType.GLN) return new Gln(atoms, id, mode,molecularSystem);
        if (type == ResidueType.ARG) return new Arg(atoms, id, mode,molecularSystem);
        if (type == ResidueType.SER) return new Ser(atoms, id, mode,molecularSystem);
        if (type == ResidueType.THR) return new Thr(atoms, id, mode,molecularSystem);
        if (type == ResidueType.VAL) return new Val(atoms, id, mode,molecularSystem);
        if (type == ResidueType.TRP) return new Trp(atoms, id, mode,molecularSystem);
        if (type == ResidueType.TYR) return new Tyr(atoms, id, mode,molecularSystem);
        return null;
    }
}
