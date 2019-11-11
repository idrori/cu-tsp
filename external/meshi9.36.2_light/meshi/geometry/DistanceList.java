package meshi.geometry;

import meshi.molecularElements.atoms.AtomCore;
import meshi.parameters.AtomType;

import java.util.ArrayList;
import java.util.Iterator;

/**
 * Created by chen on 15/12/2014.
 */
public class DistanceList extends ArrayList<Distance> {
    public final AtomCore atomOne; // The first atom of all Distances
    public final int atomOneNumber;
    public final AtomType atomOneType;
    public final boolean atomOneNowhere;
    public final boolean atomOneHeavy;


    public DistanceList(int capacity, AtomCore atomOne) {
        super(capacity);
        this.atomOne  = atomOne;
        atomOneNumber = atomOne.number;
        atomOneType   = atomOne.type();
        atomOneNowhere = atomOne.atom.nowhere();
        if (atomOne.type().isHydrogen()) atomOneHeavy = false;
        else atomOneHeavy = true;

    }
    public DistanceList(AtomCore atomOne) {
        this(DistanceLists.ROW_INITIAL_CAPACITY,atomOne);
    }


}
