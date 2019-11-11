/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import meshi.molecularElements.atoms.MolecularSystem;
import meshi.util.*;
import meshi.util.filters.Filter;

import java.util.*;

public class DistanceLists extends ArrayList<DistanceList> {
    public static final int ROW_INITIAL_CAPACITY = 10;



    public DistanceLists(int capacity) {
        super(capacity);
    }



    public void clearRows() {
         for (ArrayList row : this)
             row.clear();
    }

    public DistanceList search(int atomOneNumber) {
        return search(atomOneNumber,0,size()-1);
    }

    /*
    Binary search assuming the rows are sorted by atomNumber
     */
    private DistanceList search(int atomOneNumber,int bottom, int top) {
        int middle = (top+bottom)/2;
        DistanceList bottomRow = get(bottom);
        DistanceList middleRow = get(middle);
        DistanceList topRow    = get(top);

        if ((bottom == top) || (bottom == top-1)) {
            if (topRow.atomOneNumber == atomOneNumber) return topRow;
            if (bottomRow.atomOneNumber == atomOneNumber) return bottomRow;
            else return null;
        }
        if (middleRow.atomOneNumber == atomOneNumber) return middleRow;
        if (bottom > top) throw new RuntimeException("This is weird");
        if (middleRow.atomOneNumber > atomOneNumber) return search(atomOneNumber,bottom,middle-1);
        return search(atomOneNumber,middle+1,top);

    }

}

