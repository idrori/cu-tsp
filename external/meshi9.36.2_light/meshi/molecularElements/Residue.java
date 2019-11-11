/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.energy.beta.peptideBonds.PeptideBond;
import meshi.molecularElements.atoms.*;
import meshi.util.*;
import meshi.util.string.*;
import meshi.parameters.*;

import java.util.*;

/**
 * A generic residue.
 */
public class Residue implements ResidueCreator, Comparable, MeshiPotential, MeshiAttribute, Attributable {


    /**
     * Residue name.
     */
    public final String name;

    /**
     * A unique identifier. Replaces field number and also includes
     * this residue's chain.
     */
    private final ResidueIdentifier ID;

    /**
     * Residue number. Note that with the current setResidue of constructors,
     * it is the responsibility of the creating class to provide a meaningful number.
     **/

    /**
     * which of the 20 possible it is.
     */
    public final ResidueType type;

    public final ResidueType type() {
        return type;
    }

    /**
     * Atoms of this residue;
     */
    private AtomList atoms;

    /**
     * Co-valiant bonds in the residue.
     */
    protected AtomPairList bonds;

    /**
     * First atom of the residue (typicaly N).
     * The one conncted to the previous residue.
     */
    protected Atom head;

    /**
     * last atom of the residue (typicaly C).
     * The one conncted to the next residue.
     */
    protected Atom tail;

    protected Atom prevAtom;
    protected Atom nextAtom;
    /**
     * A list of Residue names3 (three letter code).
     */
    private static StringList names3 = new StringList();

    /**
     * A list of Residue names1 (one letter code).
     */
    private static StringList names1 = new StringList();

    /**
     * The list of residue names is done.
     */
    private static boolean namesDone = false;

    /**
     * HELIX, SHEET ect.
     * Acceptable values may be found in meshi/parameters/MeshiPotential.
     */
    private LocalStructure localStructure = new LocalStructure();
    private SecondaryStructure secondaryStructure = SecondaryStructure.UNK;
    private SecondaryStructure predictedSecondaryStructure = SecondaryStructure.UNK;
    private String accessibility = null;
    private String predictedAccessibility = null;

    private double relativeAccesibility = -1;

    public final ResidueMode mode;
    private double reliability = 0;

    private Conservation conservation = null;
    public Conservation conservation() {return conservation;}
    protected void setConservation(Conservation conservation) {this.conservation = conservation;}

    public double reliability() {
        return reliability;
    }

    public void setReliability(double reliability) {
        this.reliability = reliability;
    }

    /**
     * The PeptideBond between the residue and its presiding one.
     */
    private PeptideBond peptideBond = null;

    //----------------------------------------------------- constructors --------------------------------------------

    /**
     * <pre>
     * Use this constructor to instantiate a creator object..
     */
    public Residue(String name) {
        this(name, ResidueType.CREATOR, new ResidueIdentifier("Creator - no chain", -1, -1), null, ResidueMode.CREATOR);
    }

    public Residue(String name, ResidueType type, ResidueIdentifier id, AtomList atomList, ResidueMode mode) {
        if (atomList != null) atomList.assertSingleResidue(); //Makes sure all getAtoms belong to the same residue.
        this.name = name;
        this.ID = id;
        this.type = type;
        atoms = atomList;
        bonds = new AtomPairList();
        prevAtom = nextAtom = head = tail = null;
        if (atomList != null) atoms.setResidue(this);
        this.mode = mode;
    }

    public Residue(ResidueIdentifier id) {
        this("Dummy", ResidueType.DMY, id, null, ResidueMode.DUMMY);
    }

    public Residue(ResidueIdentifier id, String name) {
        this(name, ResidueType.DMY, id, null, ResidueMode.DUMMY);
    }

    public PeptideBond peptideBond() {
        return peptideBond;
    }

    public void setPeptideBond(PeptideBond peptideBond) {
        this.peptideBond = peptideBond;
    }
    //------------------------------------------------- methods -----------------------------------------------

    /**
     * Returns a list of the residue getAtoms.
     */
    public AtomList getAtoms() {
        return atoms;
    }

    /**
     * Returns a list of the residue bonds.
     */
    public AtomPairList bonds() {
        return bonds;
    }

    public String toString() {
        String attributes = "";
        if (secondaryStructure != SecondaryStructure.UNK) attributes += "_SS:" + secondaryStructure;
        if (accessibility != null) attributes += "_Access:" + accessibility;
        if (!ID().chain().equals(Chain.GENERIC_CHAIN_NAME))
            attributes += "_chain:" + ID().chain();
        return name + "_" + ID.number() + attributes;
    }

    public ResidueMode getMode() {
        return mode;
    }

    public String test() {
        String out = "Testing " + name + " " + comment() + "\n";

        Atom atom;
        Iterator atomIter = atoms.iterator();
        out += "getAtoms:\n";
        while ((atom = (Atom) atomIter.next()) != null)
            out += atom.comment() + "\n";

        AtomPair bond;
        Iterator bondIter = bonds.iterator();
        out += "bonds:\n";
        while ((bond = (AtomPair) bondIter.next()) != null)
            out += "\t" + bond.atom1().name + "\t" + bond.atom2().name + "\n";
        return out;
    }

    public Atom head() {
        return head;
    }

    public Atom tail() {
        return tail;
    }

    public Atom prevAtom() {
        return prevAtom;
    }

    public Atom nextAtom() {
        return nextAtom;
    }

    public void setPrevAtom(Atom atom) {
        prevAtom = atom;
    }

    public void setNextAtom(Atom atom) {
        nextAtom = atom;
    }

    public String comment() {
        return "Residue";
    }


    /**
     * residue number - 1.
     */
    public int position() {
        return ID.number() - 1;
    }

    /**
     * create a residue based isOn an getAtoms list.
     * creates a dummyResidue if the list is null or empty.
     */
    public Residue create(AtomList atoms, ResidueIdentifier id, ResidueMode mode, MolecularSystem molecularSystem) {
        if ((atoms == null) || (atoms.size() == 0)) return new Residue(id);
        Atom atom = atoms.atomAt(0);
        if (atom.molecularSystem != molecularSystem) throw new RuntimeException("This is weird.");
        String name = atom.pdbLine().residueName();
        return new Residue(name, ResidueType.type(name), id, atoms, mode);
    }


    /**
     * Compares this residue with the given residue according to
     * the order of their {@link meshi.molecularElements.ResidueIdentifier}s.
     *
     * @throws ClassCastException if the given object is not a
     *                            residue.
     * @see meshi.molecularElements.ResidueIdentifier#compareTo(Object)
     */
    public int compareTo(Object obj) {
        return ID.compareTo(((Residue) obj).ID);
    }

    /**
     * Fetches an atom by its name.
     */
    public Atom get(String name) {
        Iterator atomsIter = atoms.iterator();
        Atom atom;
        while ((atom = (Atom) atomsIter.next()) != null)
            if (atom.name.equals(name)) return atom;
        return null;
    }

    public SecondaryStructure getSecondaryStructure() {
        return secondaryStructure;
    }

    public SecondaryStructure getPredictedSecondaryStructure() {
        return predictedSecondaryStructure;
    }

    public void setSecondaryStructure(SecondaryStructure secondaryStructure) {
        this.secondaryStructure = secondaryStructure;
    }
    public void setPredictedSecondaryStructure(SecondaryStructure secondaryStructure) {
        this.predictedSecondaryStructure = secondaryStructure;
    }

    public final String getAccessibility() {
        return accessibility;
    }
    public final String getPredictedAccessibility() {
        return predictedAccessibility;
    }

    public void setAccessibility(String accessibility) {
        this.accessibility = accessibility;
    }

    public void setPredictedAccessibility(String accessibility) {
        this.predictedAccessibility = accessibility;
    }

    public Atom getAtom(String atomName) {
        if (dummy()) throw new RuntimeException("Dummy Residue - no atoms");

        for (Atom atom : atoms) {
            if (atom.name().equals(atomName))
                return atom;
        }
        return null;
    }

    public Atom ca() {
        return getAtom("CA");
    }

    public Atom cb() {
        return getAtom("CB");
    }

    public Atom cg() {
        return getAtom("CG");
    }
    public Atom cd() {
        return getAtom("CD");
    }

    public Atom carbonylC() {
        return getAtom("C");
    }

    public Atom carbonylO() {
        return getAtom("O");
    }

    public Atom amideN() {
        return getAtom("N");
    }

    public Atom amideH() {
        return getAtom("H");
    }

    public boolean dummy() {
        return type == ResidueType.DMY;
    }

    public static Atom find(AtomList atomList, BBatom bbAtom) {
        if (atomList == null) return null;
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (atom.type().bbAtom() == bbAtom) return atom;
        }
        return null;
    }

    public static Atom find(AtomList atomList, String name) {
        if (atomList == null) return null;
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (atom.name.equals(name)) return atom;
        }

        if (name.equals("HD21")) name = "1HD2";
        if (name.equals("HD22")) name = "2HD2";
        if (name.equals("HE21")) name = "1HE2";
        if (name.equals("HE22")) name = "2HE2";
        for (Iterator atoms = atomList.iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (atom.name.equals(name)) return atom;
        }

        return null;
    }


    // ------------------------ for homology modeling -----------------------------

    public void setPrevNextAtomsToNull() {
        prevAtom = null;
        nextAtom = null;

    }

    public void initiateAtoms() {
        setPrevNextAtomsToNull();
        Iterator atomsIter = atoms.iterator();
        Atom atom;
        while ((atom = (Atom) atomsIter.next()) != null) {
            atom.setResidueNumber(ID.number());
            atom.emptyBonded();
        }
    }


    public void setSecondaryStructure(char ss) {
        setSecondaryStructure(SecondaryStructure.secondaryStructure(ss));
    }

    public void setLocalStructure(LocalStructure ls){
        localStructure = new LocalStructure(ls);
    }
    public LocalStructure getLocalStructure() { return localStructure; }

    public void setAccessibility(char ac) {
        if (ac == 'A') accessibility = "ACCESSIBLE";
        if (ac == 'B') accessibility = "BURIED";
    }

    public void setRelativeccessibility(double relativeAccesibility) {
        this.relativeAccesibility = relativeAccesibility;
    }

    public double getRelativeAccesibility() {
        return relativeAccesibility;
    }
    public void defrost() {
        atoms.defrost();
    }

    /**
     * Returns the ResidueIdentifier of <code>this</code> Residue
     * or throws an exception if <code>(this != ID.RESIDUE)</code>.
     */
    public ResidueIdentifier ID() {
        return ID;
    }

    public String chain() {
        return ID().chain();
    }

    public int getChainNumber() {
         return ID().chainNumber;
    }

    public int number() {
        return ID().number();
    }

    public void setResidueInAtoms(AtomList atoms) {
        Iterator atomsIter = atoms.iterator();
        while (atomsIter.hasNext())
            ((Atom) atomsIter.next()).setResidue(this);
    }

    public boolean equals(Object other) {
        if (!(other instanceof Residue))
            return false;
        return (compareTo(other) == 0);
    }

    public int key() {
        return RESIDUE_ATTRIBUTE;
    }

    private HashMap attributes = new HashMap();

    public final void addAttribute(MeshiAttribute attribute) {
        attributes.put(attribute.key(), attribute);
    }

    public final MeshiAttribute getAttribute(int key) {
        return (MeshiAttribute) attributes.get(key);
    }

    public boolean nowhere() {
        for (Atom atom : atoms)
            if (!atom.nowhere()) return false;
        return true;
    }

    public void setNowhere() {
        for (Atom atom : atoms)
            atom.setNowhere();
    }
}
    
