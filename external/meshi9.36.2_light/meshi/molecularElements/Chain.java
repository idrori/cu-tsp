/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.molecularElements;

import meshi.geometry.Coordinates;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.MolecularSystem;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtoms;
import meshi.parameters.*;
import meshi.sequences.ResidueSequence;
import meshi.sequences.MeshiSequence;
import meshi.sequences.SequenceAlignmentCell;
import meshi.util.MeshiException;
import meshi.util.filters.Filter;

import java.util.Iterator;

/**
 * A list of residues of the same chain.
 * Note that some of the residues may be dummy.
 * Specifically, the first residue may be dummy in order to be compatible with the
 * biologists convention that the first residue is 1 (they grew up isOn FORTRAN).
 * Residue positions in the list are equal to their residue number. Holes are filed by DummyResidue objects
 */
public class Chain extends ResidueList {
    /**
     * The first Non Dummy Residue (may but also may not be number 1).
     */
    public final String name;
    public final int    number;

    public final String name() {
        return name;
    }
    public final int number() { return number;}

    public static final String GENERIC_CHAIN_NAME = " ";
    public static final int GENERIC_CHAIN_NUMBER= 0;
    public final Protein protein;

    // --------------------------------- constructors ------------------------------

    /**
     * Empty list.
     */
    public Chain(String name, int number, Protein protein) {
        super();
        this.name = name;
        this.protein = protein;
        this.number = number;
    }

    /**
     * Constructs a Chain from a sequence (a string of one letter codes).
     * The residueCreator is responsible for the interpretation of the letters to actual residues and thus
     * determines the molecular model.
     */
    public Chain(MeshiSequence sequence, ResidueCreator creator, Protein protein) {
        this(sequence, creator, GENERIC_CHAIN_NAME, GENERIC_CHAIN_NUMBER, protein);
    }

    public Chain(ResidueList residueList, Protein protein) {   //todo This constructor assumes only one chain should be modified or removed
        this(residueList.get(0).chain(), GENERIC_CHAIN_NUMBER, protein);
        int residueNumber = 0;
        for (int iResidue = 0; iResidue < residueList.size(); iResidue++) {
            Residue residue = residueList.get(iResidue);
            while (residueNumber < residue.ID().number())  {
                add(new Residue(new ResidueIdentifier(name,residueNumber, GENERIC_CHAIN_NUMBER)));
                residueNumber++;
            }
//            if ((size()> 0) && (!get(size()-1).dummy())){
//                Residue previousResidue = get(size()-1);
//                if (previousResidue.ID().number()!= residue.ID().number()-1)
//                    throw new RuntimeException("This is weird");
//                AtomPair bond = previousResidue.head().bond(residue.tail());
//                //residue.bonds().add(bond);
//            }
            add(residue);
            residueNumber = residue.ID().number+1;
        }
    }
    /**
     * Constructs a Chain from a sequence (a string of one letter codes).
     * The residueCreator is responsible for the interpretation of the letters to actual residues and thus
     * determines the molecular model.
     */
    public Chain(MeshiSequence sequence, ResidueCreator creator, String chainName, int chainNumber, Protein protein) {
        this(chainName, chainNumber, protein);
        int i;
        Residue newResidue;
        SecondaryStructure ss;
        MolecularSystem tempMolecularSystem = new MolecularSystem();
        add(new Residue(new ResidueIdentifier(chainName, 0, chainNumber))); //traditionally sequences start at 1
        if (sequence.size() < 1) throw new RuntimeException("Protein sequence length needs to be at least one");
        if (sequence.size() < 2)
            add(creator.create(getCaAtomList(sequence.getChar(0),
                    new Coordinates(),tempMolecularSystem),
                    new ResidueIdentifier(chainName, 1, chainNumber),
                    ResidueMode.SINGLE,protein.molecularSystem));
        else {
            i = 0;
            while (((i < sequence.startsIn() - 1) | (sequence.getChar(i) == 'X'))
                    & (i < sequence.size())) {
                newResidue = new Residue(new ResidueIdentifier(chainName, i + 1, chainNumber));
                add(newResidue);
                i++;
            }
            if (i == sequence.size()) throw new RuntimeException("Weired sequence");

            AtomList atomList = getCaAtomList(sequence.getChar(i - sequence.startsIn()),
                    new Coordinates(),tempMolecularSystem);
            newResidue = creator.create(atomList,
                    new ResidueIdentifier(chainName, i + 1, chainNumber),
                    ResidueMode.NTER,protein.molecularSystem);
            ss = sequence.getSs(i - sequence.startsIn() + 1);
            newResidue.setSecondaryStructure(ss);
            add(newResidue);
            i++;
           while ( i < sequence.size() + sequence.startsIn() - 1) {
                if (sequence.getChar(i - sequence.startsIn() + 1) == 'X') {
                    newResidue = new Residue(new ResidueIdentifier(chainName, i + 1, chainNumber));
                    add(newResidue);
                } else {
                    newResidue = creator.create(getCaAtomList(sequence.getChar(i - sequence.startsIn()),
                            new Coordinates(),tempMolecularSystem),
                            new ResidueIdentifier(chainName, i + 1, chainNumber),
                            ResidueMode.NORMAL,protein.molecularSystem);
                    ss = sequence.getSs(i - sequence.startsIn() + 1);
                    newResidue.setSecondaryStructure(ss);
                    add(newResidue);
                }
                i++;
            }
            newResidue = creator.create(getCaAtomList(sequence.getChar(i - sequence.startsIn()),
                    new Coordinates(),tempMolecularSystem),
                    new ResidueIdentifier(chainName, i + 1, chainNumber),
                    ResidueMode.CTER,protein.molecularSystem);
            ss = sequence.getSs(i - sequence.startsIn());
            newResidue.setSecondaryStructure(ss);
            add(newResidue);
        }

    }

    /**
     * Constructs a Chain from a list of getAtoms.
     * The residueCreator allows the use of information that is not  stored in the getAtoms themselves.
     * Note that the chain is going to havewill get a genergic name.
     */
    public Chain(AtomList atomList, ResidueCreator creator, Protein protein) {
        this(atomList, creator, GENERIC_CHAIN_NAME, protein);
    }

    /**
     * Constructs a Chain from a list of atoms.
     * The residueCreator allows the use of information that is not  stored in the atoms themselves.
     */
    public Chain(AtomList atomList, ResidueCreator creator, String chainName, Protein protein) {       //// TODO: 21/08/2017 Assumes a single chain. Should be changed or removed.
        this(chainName, GENERIC_CHAIN_NUMBER, protein);
        AtomList newAtomList = null;
        boolean first = true;
        add(new Residue(new ResidueIdentifier(chainName, 0, GENERIC_CHAIN_NUMBER))); //traditionally sequences start at 1
        if (atomList.size() == 0) throw new RuntimeException(" No Atoms in AtomList " + atomList.comment());
        for (Atom atom : atomList) {
            if (atom.residueNumber() >= 1)  {
            while (size() < atom.residueNumber()) {
                if (newAtomList == null) add(new Residue(new ResidueIdentifier(chainName, size(), GENERIC_CHAIN_NUMBER)));
                else {
                    if (first) {
                        boolean success = createAndAddResidue(creator, newAtomList, ResidueMode.NTER);
                        if (success)first = false;
                    } else {
                        createAndAddResidue(creator, newAtomList, ResidueMode.NORMAL);
                    }
                    newAtomList = null;
                }
            }
            if (size() == atom.residueNumber()) {
                if (newAtomList == null) newAtomList = new AtomList(atom.molecularSystem);
                    newAtomList.add(atom);
                }
            }
        }
        createAndAddResidue(creator, newAtomList, ResidueMode.CTER);
        sort();
    }


//    public Chain(AtomList atomList, ResidueCreator creator, String chainName, Protein protein, MolecularSystem current, MolecularSystem newms) {
//        this(chainName, protein);
//        AtomList newAtomList = null;
//        boolean first = true;
//        add(new Residue(new ResidueIdentifier(chainName, 0))); //traditionally sequences start at 1
//        if (atomList.size() == 0) throw new RuntimeException(" No Atoms in AtomList " + atomList.comment());
//        for (Atom atom : atomList) {
//            while (size() < atom.residueNumber()) {
//                if (newAtomList == null) add(new Residue(new ResidueIdentifier(chainName, size())));
//                else {
//                    if (first) {
//                        createAndAddResidue(creator, newAtomList, ResidueMode.NTER);
//                        first = false;
//                    } else {
//                        createAndAddResidue(creator, newAtomList, ResidueMode.NORMAL);
//                    }
//                    newAtomList = null;
//                }
//            }
//            if (size() == atom.residueNumber()) {
//                if (newAtomList == null) newAtomList = new AtomList(atom.molecularSystem);
//                newAtomList.add(atom);
//            }
//        }
//        createAndAddResidue(creator, newAtomList, ResidueMode.CTER);
//        sort();
//    }

    private static AtomList getCaAtomList(char c, Coordinates coordinates,MolecularSystem molecularSystem) {
        AtomType type = ResidueType.type(c).caType();
        if (type == AtomType.XXX) throw new RuntimeException("Problem in getCaAtomList" +
                " cannot convert character \"" + c + "\" to a residue type.");
        Atom ca = new Atom("CA", null, type, coordinates, 1, new Double(0),molecularSystem);
        AtomList out = new AtomList(molecularSystem);
        out.add(ca);
        return out;
    }

    private boolean createAndAddResidue(ResidueCreator creator, AtomList newAtomList, ResidueMode mode) {
        Atom ca = ResidueExtendedAtoms.find(newAtomList, BBatom.CA);
        if (ca == null) return false;
        ResidueIdentifier ri = new ResidueIdentifier(name(), size(), number);
        Residue residue;
        residue = creator.create(newAtomList, ri, mode,protein.molecularSystem);
        add(residue);
        residue.setResidueInAtoms(newAtomList);
        return true;
    }

/*
	private void createAndAddResidue(ResidueCreator creator, int mode) {
		Residue residue = creator.create(newAtomList,size(),Residues.NORMAL);
	}
*/


    //------------------------------------------- methods -----------------------------------

    /**
     * Fetch a residue by its position in the list.
     */
    public Residue residueAt(int index) {
        return get(index);
    }

    /**
     * Fetches a residue by its residue identifier.
     */
    public Residue residue(ResidueIdentifier residueID) {
        Residue residue;
        Iterator residueIterator = iterator();
        while (residueIterator.hasNext()) {
            residue = (Residue) residueIterator.next();
            if (residue.ID().equals(residueID))
                return residue;
        }
        return null;
    }

    /**
     * For users who need this signature.
     */
// 	public Residue residue(int residueNumber) {
// 		return residue(new Residue(residueNumber).ID());
// 	}

// Previous implementation:
/*
     public Residue residue(int residueNumber) {
         if (residueNumber >= size()) return null;
	 Residue out = get(residueNumber);
	 if (residueNumber < 0) 
	     throw new RuntimeException("Weird residueNumber "+residueNumber);
	 if (out.number != residueNumber) 
	     throw new RuntimeException("Resdiue "+out+" in position "+residueNumber+"\n"+
					"in "+this);
	 return out;
	 }
*/

// 	Residue key = new Residue(residueNumber);
// 	int index;
// 	try {
// 	    index = binarySearch(key);
// 	}
// 	catch (Exception ex) {
// 	    System.out.println("************ ERROR **************");
// 	    print(10);
// 	    ex.printStackTrace();
// 	    throw new RuntimeException("A problem with residueNumber "+residueNumber+"\n"+
// 				       "key = "+key+
// 				       "size = "+size()+"\n"+
// 				       "residues(residueNumber) = "+get(residueNumber)+"\n"+
// 				       ex+"\n");
// 	}
// 	if (index < 0) return null;
// 	return residueAt(index);

    /**
     * Extracts the residues accepted by the filter.
     */
    public ResidueList filter(Filter residueFilter) {
        ResidueList newList = new ResidueList();
        for (Residue residue : this) {
            if (residueFilter.accept(residue))
                newList.add(residue);
        }
        return newList;
    }

    static class IsResidue implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof Residue);
        }
    }


    public int numberOfNonDummyResidues() {
        return getNonDummyResidues().size();
    }

    public ResidueList getNonDummyResidues() {
        ResidueList out = new ResidueList();
        for (Residue residue : this)
            if (!residue.dummy()) out.add(residue);
        return out;
    }


    public String toString() {
        String out = "Chain " + name + ":\n";
        Iterator residues = iterator();
        Residue residue = (Residue) residues.next(); // remove the first dummy residue
        while (residues.hasNext()) {
            residue = (Residue) residues.next();
            if (residue.dummy()) out += SequenceAlignmentCell.WILDCARD_CHAR;
            else out += residue.type().nameOneLetter();
        }
        return out;
    }

    public ResidueSequence sequence() {
        return new ResidueSequence(this, protein.name() + "_" + name());
    }
    public String sequenceWithXs() {
        String out = "";
        boolean first = true;
        for (Residue residue : this) {
            if (residue.dummy())
                if (!first)
                    out +="X";
                else
                    first = false;
            else
                out += residue.type().nameOneLetter();
        }
        return out;
    }
    public AtomList atoms() {
        return new AtomList(this,protein.molecularSystem);
    }

//     public void setChain(String chain) {
// //		assertSingleChain();
// 		Iterator residueIter = iterator();
// 		Residue residue;
// 		while (residueIter.hasNext()) {
// 			residue = (Residue) residueIter.next();
// 			residue.setChain(chain);
// 		}
// 	}

    /**
     * Asserts that <code>this</code> chain is comprised of
     * residues with the same chain string.
     *
     * @throws meshi.util.MeshiException if the assertion fails.
     */
    public void assertChain() {
        for (Iterator residueIter = iterator(); residueIter.hasNext();) {
            Residue residue = (Residue) residueIter.next();
            if (!residue.ID().chain().equals(name()))
                throw new MeshiException("Chain.assertChain -- " +
                        "It appears that this Chain " + this + " contains residue " + residue + "  of other chain.");
        }
    }

    public ResidueList missingResidues() {
        ResidueList out = new ResidueList();
        boolean nonDummyFound = false;
        for (Iterator residues = iterator(); residues.hasNext() & (!nonDummyFound);) {
            Residue residue = (Residue) residues.next();
            if (!residue.dummy()) nonDummyFound = true;
        }
        if (nonDummyFound) {
            for (Iterator residues = iterator(); residues.hasNext() & (!nonDummyFound);) {
                Residue residue = (Residue) residues.next();
                if (residue.dummy()) out.add(residue);
            }
        }
        return out;
    }

    public AtomList nowhereAtoms(MolecularSystem molecularSystem) {
        AtomList out = new AtomList(molecularSystem);
        for (Iterator atoms = atoms().iterator(); atoms.hasNext();) {
            Atom atom = (Atom) atoms.next();
            if (atom.nowhere()) out.add(atom);
        }
        return out;
    }

    public Residue firstNonDummyResidue() {
        int index = 0;
        for (Iterator residues = iterator(); residues.hasNext();) {
            Residue residue = (Residue) residues.next();
            if (!residue.dummy()) {
                if (residue.number() != index) throw new RuntimeException("This is weird");
                return residue;
            }
            index++;
        }
        return null;
    }

}

