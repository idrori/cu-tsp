/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.alignmentTools;

import java.util.*;

import meshi.util.*;
import meshi.util.overlap.*;
import meshi.util.filters.*;
import meshi.molecularElements.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.*;
import meshi.geometry.*;
import meshi.geometry.rotamers.*;
import meshi.molecularElements.extendedAtoms.*;


public class AlignmentTools implements KeyWords {
    public static void thread(Protein prot, int begins, int ends, int move) {
        AtomList a1, a2;
        Atom at1, at2;
        for (int res = begins; res <= ends; res++) {
            a1 = prot.residue(res).getAtoms().filter(new AtomList.BackboneFilter());
            a2 = prot.residue(res + move).getAtoms().filter(new AtomList.BackboneFilter());
            for (int c = 0; c < a1.size(); c++) {
                at1 = a1.atomAt(c);
                at2 = getAtom(a2, at1.name());
                if (at2 != null)
                    at1.setXYZ(at2.x(), at2.y(), at2.z());
            }
        }
    }

    public static void outOfTheWay(Protein prot, int resnum, double shift) {
        double x, y, z, cmx = 0, cmy = 0, cmz = 0; // center of mass x, y and z
        Atom atom;
        Iterator iter = prot.atoms().iterator();
        while ((atom = (Atom) iter.next()) != null) {
            cmx += atom.x();
            cmy += atom.y();
            cmz += atom.z();
        }
        cmx /= prot.atoms().size();
        cmy /= prot.atoms().size();
        cmz /= prot.atoms().size();

        AtomList a1 = prot.residue(resnum).getAtoms();
        for (int c = 0; c < a1.size(); c++) {
            x = a1.atomAt(c).x();
            y = a1.atomAt(c).y();
            z = a1.atomAt(c).z();
            double norm = Math.sqrt((x - cmx) * (x - cmx) +
                    (y - cmy) * (y - cmy) +
                    (z - cmz) * (z - cmz));
            a1.atomAt(c).setXYZ((cmx + (norm + shift) / norm * (x - cmx)),
                    (cmy + (norm + shift) / norm * (y - cmy)),
                    (cmz + (norm + shift) / norm * (z - cmz)));
        }
    }

    public static Atom getAtom(AtomList al, String name) {
        for (int c = 0; c < al.size(); c++)
            if (al.atomAt(c).name().equals(name))
                return al.atomAt(c);
        return null;
    }

    public static Atom getAtom(AtomList al, int resNum, String name) {
        for (int c = 0; c < al.size(); c++)
            if ((al.atomAt(c).residueNumber() == resNum) && al.atomAt(c).name().equals(name))
                return al.atomAt(c);
        return null;
    }
}




