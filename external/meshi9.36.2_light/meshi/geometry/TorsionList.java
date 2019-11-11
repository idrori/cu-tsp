/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.util.*;
import meshi.util.filters.*;
import meshi.molecularElements.atoms.*;
import meshi.molecularElements.*;
import meshi.energy.simpleEnergyTerms.GoodBonds;

import java.util.*;

/**
 **/
public class TorsionList extends ArrayList<Torsion> implements Updateable {
    private int numberOfUpdates = 0;


    /**
     * An Torsion list based isOn a angle list
     */
    public TorsionList() {
        super();
    }

    public TorsionList(AngleList angles, DistanceMatrix distanceMatrix) {
        super();
        Angle angle1, angle2;
        Torsion tor;
        for (int iAngle1 = 0; iAngle1 < angles.size(); iAngle1++) {
            angle1 = angles.get(iAngle1);
            for (int iAngle2 = 0; iAngle2 < iAngle1; iAngle2++) {
                angle2 = angles.get(iAngle2);
                if ((angle1.sharedAtomPair(angle2) != null) &&      // There is a shared pair and the angles are not the same but in opposite directions
                        !((angle1.atom1 == angle2.atom3) && (angle1.atom2 == angle2.atom2) &&
                                (angle1.atom3 == angle2.atom1))) {
                    tor = getTorsion(angle1, angle2, distanceMatrix);
                    if (isNamed(tor))
                        add(tor);
                    tor = getTorsion(angle2, angle1, distanceMatrix);
                    if (isNamed(tor))
                        add(tor);
                }
            }
        }
    }

    public Torsion getTorsion(Angle angle1, Angle angle2, DistanceMatrix distanceMatrix) {
//        return new Torsion(angle1, angle2, distanceMatrix);
        return new QuickAndDirtyTorsion(angle1, angle2, distanceMatrix);
    }

    public void update(int numberOfUpdates) throws UpdateableException {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            int size = size();
            for (int i = 0; i < size; i++) {
                torsionAt(i).update(numberOfUpdates);
            }
            this.numberOfUpdates++;
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with TorsionList.update(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }


    public Torsion torsionAt(int i) {
        return (Torsion) get(i);
    }

    public AtomList atomList() {
        if (size() < 1) throw new RuntimeException("Cannot create an atom listf from an empty list of torsions");
        AtomList list = new AtomList(get(0).atom1.molecularSystem);
        Atom atom1, atom2, atom3, atom4;
        for (Torsion torsion : this) {
            atom1 = torsion.atom1;
            atom2 = torsion.atom2;
            atom3 = torsion.atom3;
            atom4 = torsion.atom4;
            if (!list.contains(atom1)) {
                list.add(atom1);
            }
            if (!list.contains(atom2)) {
                list.add(atom2);
            }
            if (!list.contains(atom3)) {
                list.add(atom3);
            }
            if (!list.contains(atom4)) {
                list.add(atom4);
            }
        }
        return list;
    }

    public boolean equivalentExists(Torsion findMe) {
        Iterator torsions = iterator();
        for (Torsion torsion : this) {
            if (torsion.equivalent(findMe)) return true;
        }
        return false;
    }

    public TorsionList filterEquivalents() {
        TorsionList out = new TorsionList();
        for (Torsion torsion : this) {
            if (isPhi(torsion)) out.add(torsion);
            else if (isPsi(torsion)) out.add(torsion);
            else if (isOmega(torsion)) out.add(torsion);
            else if (isChi1(torsion)) out.add(torsion);
            else if (isChi2(torsion)) out.add(torsion);
            else if (isChi3(torsion)) out.add(torsion);
            else if (isChi4(torsion)) out.add(torsion);
            else if (isNimp(torsion)) out.add(torsion);
            else if (isCimp(torsion)) out.add(torsion);
        }
        for (Torsion torsion : this) {
            if (!out.equivalentExists(torsion)) out.add(torsion);
        }
        return out;
    }

    /**
     * Returns a sub-list that is accepted by the parameter
     */
    public TorsionList chi1Filter() {
        TorsionList out = new TorsionList();
        for (Torsion torsion : this) {
            if (isChi1(torsion)) out.add(torsion);
        }
        return out;
    }

    /**
     * Returns a sub-list that is has a known name
     */
    public TorsionList namedFilter() {
        TorsionList out = new TorsionList();
        for (Torsion torsion : this) {
            if (isNamed(torsion)) out.add(torsion);
        }
        return out;
    }

    public static boolean isNamed(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("") != 0)
            return true;
        else
            return false;
    }

    public static boolean isPhi(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("PHI") == 0)
            return true;
        else
            return false;
    }

    public static boolean isPsi(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("PSI") == 0)
            return true;
        else
            return false;
    }

    public static boolean isOmega(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("OMG") == 0)
            return true;
        else
            return false;
    }

    public static boolean isCimp(Torsion torsion) {
        if (torsion.atom1.name.equals("CA") &
                torsion.atom2.name.equals("N") &
                torsion.atom3.name.equals("C") &
                torsion.atom4.name.equals("O")) return true;
        return false;
    }

    public static boolean isNimp(Torsion torsion) {
        if (torsion.atom1.name.equals("H") &
                torsion.atom2.name.equals("C") &
                torsion.atom3.name.equals("N") &
                torsion.atom4.name.equals("CA")) return true;
        return false;
    }

    public static boolean isChi1(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("CHI1") == 0)
            return true;
        else
            return false;
    }

    public static boolean isChi2(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("CHI2") == 0)
            return true;
        else
            return false;
    }

    public static boolean isChi3(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("CHI3") == 0)
            return true;
        else
            return false;
    }

    public static boolean isChi4(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("CHI4") == 0)
            return true;
        else
            return false;
    }

    public static boolean isChi5(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("CHI5") == 0)
            return true;
        else
            return false;
    }

    public static boolean isOOP(Torsion torsion) {
        if (torsion.getTorsionName().compareTo("OOP") == 0)
            return true;
        else
            return false;
    }


    public static class FilterOOP implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isOOP(tor);
        }
    }

    public static class FilterPhi implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isPhi(tor);
        }
    }

    public static class FilterPsi implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isPsi(tor);
        }
    }

    public static class FilterChi1 implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isChi1(tor);
        }
    }

    public static class FilterChi2 implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isChi2(tor);
        }
    }

    public static class FilterChi3 implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isChi3(tor);
        }
    }

    public static class FilterChi4 implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return isChi4(tor);
        }
    }

    public static class FilterSideChain implements Filter {
        public boolean accept(Object obj) {
            Torsion tor = (Torsion) obj;
            return (isChi1(tor) || isChi2(tor) || isChi3(tor) || isChi4(tor));
        }
    }


    public static TorsionList createTorsionList(Protein protein, DistanceMatrix distanceMatrix) {
        AtomPairList bondList = (AtomPairList) protein.bonds().filter(new GoodBonds());
        AngleList angleList = new AngleList(bondList, distanceMatrix);
        return new TorsionList(angleList, distanceMatrix);
    }

    public static TorsionList createQuickAndDirtyTorsionList(Protein protein, DistanceMatrix distanceMatrix) {
        AtomPairList bondList = (AtomPairList) protein.bonds().filter(new GoodBonds());
        AngleList angleList = new AngleList(bondList, distanceMatrix);
        return new QuickAndDirtyTorsionList(angleList, distanceMatrix);
    }

    public void freeze() {
        for (Torsion torsion : this) {
            torsion.freeze();
        }
    }

    public TorsionList filter(Filter filter) {
        TorsionList out = new TorsionList();
        for (Torsion t : this) {
            if (filter.accept(t)) out.add(t);
        }
        return out;
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
        for (Torsion torsion : this)
            torsion.setNumberOfUpdates(numberOfUpdates);
    }

}

