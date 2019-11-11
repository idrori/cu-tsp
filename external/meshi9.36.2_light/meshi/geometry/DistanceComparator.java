/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry;

import java.util.Comparator;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 18/01/2010
 * Time: 10:47:17
 * To change this template use File | Settings | File Templates.
 */
public class DistanceComparator implements Comparator {
    public static final DistanceComparator distanceComparator = new DistanceComparator();
    private Distance d1, d2;

    public int compare(Object o1, Object o2) {
        d1 = (Distance) o1;
        d2 = (Distance) o2;
        try {
            if (d1.distance() < d2.distance()) return -1;
        }
        catch (RuntimeException ex) {
            System.out.println("d1 = " + d1 + "       d2 = " + d2);
            throw ex;
        }
        if (d1.distance() == d2.distance()) return 0;
        return 1;
    }
}
