package meshi.energy.KB2013;

import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.AtomPair;

import java.util.Arrays;
import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 30/07/13
 * Time: 10:23
 * To change this template use File | Settings | File Templates.
 */
public class AtomPairIterator implements Iterator<AtomPair> {
    private Atom[] atoms;
    private int iAtom,jAtom,nAtoms,bondedDepth;
    private AtomPair[] current = new AtomPair[2];
    private AtomList bonded;
    private Atom atomI,atomJ;
    boolean hasNext;
    int i = 0;
    public AtomPairIterator(AtomList atoms, int bondedDepth) {
        this.atoms       = atoms.toArray(new Atom[1]);
        Arrays.sort(this.atoms);
        iAtom            = 0;
        jAtom            = 0;
        nAtoms           = atoms.size();
        this.bondedDepth = bondedDepth;
        bonded           = new AtomList(atoms.molecularSystem);
        iAtom = skipHydrogenAndNowhereAtoms(iAtom);
        atomI = this.atoms[iAtom];
        current[0] = new AtomPair(atomI.molecularSystem);
        DistanceMatrix.getBonded(atomI, bondedDepth, bonded, atomI.number());
        hasNext = getNext(current);
    }

    private int skipHydrogenAndNowhereAtoms(int start) {
        while((start<nAtoms) && (atoms[start].isHydrogen() || atoms[start].nowhere()))start++;
        return start;
    }

    public boolean hasNext() {return hasNext;}
    public AtomPair next() {
        current[1] = current[0];
        hasNext      = getNext(current);
        return current[1];
    }
    public void remove() {
        throw new RuntimeException(" Remove is not implemented.");
    }

    private boolean getNext(AtomPair[] current) {
        if (iAtom >= nAtoms) return false;
        if((jAtom < nAtoms)&&((jAtom <= iAtom) || (atoms[jAtom].residueNumber()<=atoms[iAtom].residueNumber()))){
            jAtom++;
            jAtom = skipHydrogenAndNowhereAtoms(jAtom);
            if (jAtom >= nAtoms)return false;
            return getNext(current);
        }
        if (jAtom == nAtoms) {
            iAtom++;
            if (iAtom == nAtoms) return false;
            iAtom = skipHydrogenAndNowhereAtoms(iAtom);
            if (iAtom == nAtoms) return false;
            atomI = atoms[iAtom];
            bonded.clear();
            DistanceMatrix.getBonded(atomI, bondedDepth, bonded, atomI.number());
            jAtom = iAtom+1;
            return getNext(current);
        }
        jAtom = skipHydrogenAndNowhereAtoms(jAtom);
        if (jAtom == nAtoms) return getNext(current);
        atomJ = atoms[jAtom];
        jAtom++;
        if (jAtom == nAtoms) return getNext(current);
        if(bonded.contains(atomJ)) return getNext(current);
        current[0].set(atomI,atomJ);
        return true;
    }
}
