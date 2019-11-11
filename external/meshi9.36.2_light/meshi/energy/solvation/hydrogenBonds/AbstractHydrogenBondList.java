/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvation.hydrogenBonds;

import java.util.Vector;

import meshi.geometry.*;
import meshi.molecularElements.atoms.*;
import meshi.parameters.AtomType;
import meshi.util.Updateable;
import meshi.util.UpdateableException;
import meshi.util.Utils;


/**
 * A super class with the most general treatment of a list of hydrogen bonds in an atom list.
 * A specific implementation should put content to these abstract methods according to the specific
 * implementaion of the hydrogen bonds. The constructor gets a distance matrix object, and an
 * AtomList object of all the getAtoms where bonds are searched. The class is centered around the
 * "bondList" vector that stores all the hydrogen bond objects relevent to the atom list. Whenever
 * a new hydrogen bond appears in the non-bonded-list, it is added to this vector. Hydrogen bonds that are
 * no longer in the non-bonded list are occasionaly removed from this vector. The parameter
 * 'refreshVectorEvery' in the constructor determine how often this 'clean-up' is done. The bondList
 * vector is non-redundent.
 * <p/>
 * <p/>
 * IPORTANT IMPORTANT IMPORTANT
 * NOTE: The distance matrix must be in an updated state before any of this class's
 * methods are called. If it is not up to date, the results would be meaningless!!!
 * This class does not update the distance matrix isOn its own at any time!!
 */
public abstract class AbstractHydrogenBondList implements Updateable {

    // Fields:
    // -------
    protected int maxAtomNum;
    protected int atomListSize;
    protected DistanceMatrix dm = null;
    protected AtomList atomList = null;
    private DistanceLists newToNonBonded; // A pointer to this list is added to the DistanceMatrix's list of lists (see line 74). This way, only new distances are recognized.
    protected Vector<AbstractHydrogenBond> bondList;
    protected int[] lut;     // The look-up table (lut) converts the atom internal number (field of Atom) to its index in the atom list given to the constructor.
    protected AbstractHydrogenBond[][] lutHB; // lutHB[i] is an array of pointers to all the hydrogen bond objects where atom number i is either their donor or acceptor. It is null for non-polar getAtoms, or getAtoms not participating in a hydrogen bond.
    private int numberOfUpdates = 0; // for the Updateable interface
    private int refreshVectorEvery = 200;


    public AbstractHydrogenBondList() {
    }

    public AbstractHydrogenBondList(DistanceMatrix dm, AtomList atomList, int refreshVectorEvery) {
        this.dm = dm;
        this.atomList = atomList;
        this.refreshVectorEvery = refreshVectorEvery;
        bondList = new Vector<AbstractHydrogenBond>();

        // Building the look-up tables
        maxAtomNum = -1;
        atomListSize = atomList.size();
        for (int c = 0; c < atomListSize; c++) {
            if (atomList.atomAt(c).number() > maxAtomNum)
                maxAtomNum = atomList.atomAt(c).number();
        }
        lut = new int[maxAtomNum + 1];
        lutHB = new AbstractHydrogenBond[maxAtomNum + 1][];
        for (int c = 0; c < maxAtomNum; c++) {
            lut[c] = -1;
        }
        for (int c = 0; c < atomListSize; c++) {
            lut[atomList.atomAt(c).number()] = c;
        }

        // This part is needed so that we can get the new distances that were added
        // to the non-bonded list in the last update.
        newToNonBonded = new DistanceLists(maxAtomNum);
        dm.energyTermsDistanceLists().add(newToNonBonded);

        buildSpecificStructures();

    }


    public void update(int updateNumber) throws UpdateableException {
        if (updateNumber == numberOfUpdates + 1) {
            if ((updateNumber % refreshVectorEvery == 0) && (updateNumber > 0))
                update(true, updateNumber);
            else
                update(false, updateNumber);
            this.numberOfUpdates++;
        } else if (updateNumber != this.numberOfUpdates)
            throw new RuntimeException("Something weird with the update in AbstractHydrogenBondList.java\n" +
                    "The update number is " + updateNumber + " while the numberOfUpdates in the class is " + numberOfUpdates);
    }


    /**
     * The 'toRefreshVector' parameter determines if broken H-bonds are removed from the vector.
     */
    protected void update(boolean toRefreshVector, int updateNumber) throws UpdateableException {
        addNewBondsToList();
        for (AbstractHydrogenBond hb : bondList)
            hb.updateHBvalueAndDerivatives();
        if (toRefreshVector)
            removeBrokenHB();
    }


    /**
     * This methods add new h-bonds to the bondList vector. Care is taken to keep the vector non-redundent.
     */
    private void addNewBondsToList() {
        DistanceLists newNonBondedList;
        newNonBondedList = dm.newNonBondedList();
        AbstractHydrogenBond hb;
        for (DistanceList newNonBondedRow : newNonBondedList) {
            Atom atom1 = newNonBondedRow.atomOne.atom;
            if (atom1.type().isCarbon() || atom1.type().isHydrogen() ) continue; // We seek oxygen and nitrogen getAtoms.
            for (Distance distance: newNonBondedRow)  {
                Atom atom2 = distance.atom2();
                if ( atom2.type().isCarbon() ||  atom2.type().isHydrogen()) continue;
                // Possible HB Ahoy!! (atom1 and atom2 are capable of forming a hydrogen bond)
                // We first need to check that this pair of polars is not already in the vector.
                if (findBondByPolars(atom1, atom2) == null) {
                    if (distance.mode() == DistanceMode.INFINITE)
                        throw new RuntimeException("Weird distance in NBlist1 " + distance + " " + distance.atom1 + " " + distance.atom2);
                    hb = createHBfromPolars(atom1, atom2);
                    if (hb != null) {
                        bondList.add(hb);
                        addToLUTvec(lutHB, atom1, hb);
                        addToLUTvec(lutHB, atom2, hb);
                    }
                }
            }
        }
        newToNonBonded.clearRows();
    }


    /**
     * This method remove from the vector Hydrogen Bond objects, whose distance between polars is
     * so large, they no longer appears in the non-bonded list of the distance matrix.
     */
    private void removeBrokenHB() {
        Vector<AbstractHydrogenBond> tmpVector = new Vector<AbstractHydrogenBond>();
        AbstractHydrogenBond[][] tmpLutHB = new AbstractHydrogenBond[lutHB.length][];
        AbstractHydrogenBond hb;
        DistanceLists distanceLists = dm.nonBondedList();
        for (DistanceList row : distanceLists) {
            Atom atom1 = row.atomOne.atom;
            if (atom1.type().isCarbon() || atom1.type().isHydrogen()) continue;
            for (Distance distance :row) {
                Atom atom2 = distance.atom2();
                if (atom1.type().isCarbon() || atom1.type().isHydrogen()) continue;
            // Possible HB Ahoy!! (atom1 and atom2 are capable of forming a hydrogen bond)
                hb = findBondByPolars(atom1, atom2);
                if (hb != null) {
                    tmpVector.add(hb);
                    addToLUTvec(tmpLutHB, atom1, hb);
                    addToLUTvec(tmpLutHB, atom2, hb);
                }
            }
        }
        bondList.clear();
        bondList = tmpVector;
        lutHB = tmpLutHB;
    }

    /**
     * This method returns a pointer to the hydrogen bond that exists between two POLAR getAtoms (hydrogens are not consider polar)
     * If no such object exists then null is returned.
     */
    public AbstractHydrogenBond findBondByPolars(Atom atom1, Atom atom2) {
        return findBondByPolars(atom1.number(), atom2.number());
    }

    public AbstractHydrogenBond findBondByPolars(int atom1number, int atom2number) {
        if (lutHB[atom1number] != null) {
            int length = lutHB[atom1number].length;
            for (int c = 0; c < length; c++)
                if ((lutHB[atom1number][c].getFirstPolar().number() == atom2number) ||
                        (lutHB[atom1number][c].getSecondPolar().number() == atom2number))
                    return lutHB[atom1number][c];
        }
        return null;
    }


    /**
     * Auxilarry method: Updating a vector, 'vec' similar to the lutHB field, with the new HB
     * The reason we do not access lutHB directly, is that this method is also called from the
     * 'removeBrokenHB' method, and there a different vector need updating.
     */
    private void addToLUTvec(AbstractHydrogenBond[][] vec, Atom atom, AbstractHydrogenBond hb) {
        if (vec[atom.number()] == null) {
            vec[atom.number()] = new AbstractHydrogenBond[1];
            vec[atom.number()][0] = hb;
        } else {
            AbstractHydrogenBond[] tmp = new AbstractHydrogenBond[vec[atom.number()].length + 1];
            for (int c = 0; c < vec[atom.number()].length; c++)
                tmp[c] = vec[atom.number()][c];
            tmp[vec[atom.number()].length] = hb;
            vec[atom.number()] = tmp;
        }
    }

    public void print() {
        for (int c = 0; c < bondList.size(); c++) {
            System.out.println(bondList.get(c));
        }
    }


    /**
     * This method should build any implementation-specific data structures needed.
     */
    protected abstract void buildSpecificStructures();

    /**
     * This method should build an implementaion specific HB from two polar getAtoms given as parameters.
     * If the creation is not possible then null is returned.
     */
    protected abstract AbstractHydrogenBond createHBfromPolars(Atom atom1, Atom atom2);

    public final AbstractHydrogenBond[][] lutHB() {
        return lutHB;
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }
}


  

