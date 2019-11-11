/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.geometry.putH;

import meshi.molecularElements.atoms.*;

import java.util.*;

public class PutHposLog {
    public final Atom hydrogen, heavy;
    public final int numberOflegandsOfHeavy, numberOfLegandsWithCoordinatesOfHeavy;

    public PutHposLog(Atom hydrogen, Atom heavy,
                      int numberOflegandsOfHeavy, int numberOfLegandsWithCoordinatesOfHeavy) {
        this.hydrogen = hydrogen;
        this.heavy = heavy;
        this.numberOflegandsOfHeavy = numberOflegandsOfHeavy;
        this.numberOfLegandsWithCoordinatesOfHeavy = numberOfLegandsWithCoordinatesOfHeavy;
    }
}
