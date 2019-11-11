/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements.ca;

import meshi.geometry.Coordinates;
import meshi.geometry.Distance;
import meshi.geometry.FreeDistance;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueIdentifier;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.parameters.AtomType;
import meshi.parameters.BBatom;
import meshi.parameters.ResidueMode;
import meshi.parameters.ResidueType;

import java.util.Iterator;


public class CaResidue extends Residue {
    public static final CaResidue creator = new CaResidue("creator");
    public final Atom CA;


    /**
     * <pre>
     * Use this constructor to instantiate a creator object..
     */
    public CaResidue(String name) {
        super(name);
        CA = null;
    }

    public CaResidue(String name, ResidueType type, ResidueIdentifier id, AtomList atomList, ResidueMode mode,MolecularSystem molecularSystem) {
        super(name, type, id, getNewAtomList(atomList, type,molecularSystem), mode);
        CA = getAtoms().get(0);
        CA.setResidue(this);
        head = tail = CA;
    }

    public static AtomList getNewAtomList(AtomList oldAL, ResidueType type,MolecularSystem molecularSystem) {
        Atom old = find(oldAL, BBatom.CA);
        if (old == null) throw new RuntimeException("CA undefined");
        return new AtomList(new Atom("CA", null, type.caType(), new Coordinates(old), old.occupancy(), old.temperatureFactor(),molecularSystem));
    }


    public int getContacts(ResidueList residueList, double contactDistance) {
        Iterator residues = residueList.iterator();
        Residue residue;
        double dis;
        int contacts = 0;
        while ((residue = (Residue) residues.next()) != null) {
            if (!(residue.dummy())) {
                dis = (distance((CaResidue) residue)).distance();
                if (dis < contactDistance) contacts++;
            }
        }
        return contacts;
    }


    public void setCoordinates(CaResidue from) {
        moveTo(from.CA.x(),
                from.CA.y(),
                from.CA.z());
    }

    public void setCoordinates(Atom from) {
        moveTo(from.x(),
                from.y(),
                from.z());
    }

    public Distance distance(CaResidue other) {
        return new FreeDistance(CA, other.CA);
    }


    public void moveTo(double x, double y, double z) {
        CA.setXYZ(x, y, z);
    }

    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode,MolecularSystem molecularSystem) {
        if ((atoms == null) || (atoms.size() == 0)) return new Residue(id);
        Atom atom = find(atoms, BBatom.CA);
        AtomType atomType = atom.type();
        ResidueType residueType = ResidueType.type(atomType);
        String residueName = residueType.nameThreeLetters();
        return new CaResidue(residueName, residueType, id, atoms, mode,molecularSystem);
    }
}
