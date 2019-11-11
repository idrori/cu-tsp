/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.beta.peptideBonds;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 24/01/2010
 * Time: 09:05:41
 * To change this template use File | Settings | File Templates.
 */
public class PeptideBondDistanceList extends ArrayList<PeptideBondDistance> {
    public PeptideBondDistanceList(int capacity) {
        super(capacity);
    }

    public static PeptideBondDistanceList sort(PeptideBondDistanceList inList) {
        PeptideBondDistanceList outList = new PeptideBondDistanceList(inList.size());
        Object[] array = inList.toArray();
        Arrays.sort(array, PeptideBondDistanceComparator.comparator);
        for (Object o : array)
            outList.add((PeptideBondDistance) o);
        //for (PeptideBondDistance pbd : outList) System.out.println("Sorted "+pbd);
        return outList;
    }

    private static class PeptideBondDistanceComparator implements Comparator {
        public static final PeptideBondDistanceComparator comparator = new PeptideBondDistanceComparator();

        public int compare(Object o1, Object o2) {
            PeptideBondDistance pbd1 = (PeptideBondDistance) o1;
            PeptideBondDistance pbd2 = (PeptideBondDistance) o2;
            if (pbd1.distance() < pbd2.distance()) return -1;
            if (pbd1.distance() > pbd2.distance()) return 1;
            return 0;
        }
    }

    public PeptideBondDistance min() {
        PeptideBondDistance out = null;
        for (PeptideBondDistance pbd : this) {
            if (out == null) out = pbd;
            else if (out.distance() > pbd.distance()) out = pbd;
        }
        return out;
    }
}
