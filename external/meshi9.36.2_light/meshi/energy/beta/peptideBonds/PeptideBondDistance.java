/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.beta.peptideBonds;

import meshi.energy.beta.peptideBonds.PeptideBond;
import meshi.geometry.FreeDistance;


/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 22/01/2010
 * Time: 17:41:40
 * To change this template use File | Settings | File Templates.
 */
public class PeptideBondDistance {
    public static final double ALPHA = -5;
    public final PeptideBond pb1, pb2;
    public final FreeDistance pb1Hpb2O, pb2Hpb1O, hhDistance, ooDistance;
    private double d12;
    private double d21;
    private double distance, dDistanceDx12, dDistanceDy12, dDistanceDz12, dDistanceDx21, dDistanceDy21, dDistanceDz21;
    double d12ALPHA, expD12ALPHA, dExpD12ALPHAdD12;
    double d21ALPHA, expD21ALPHA, dExpD21ALPHAdD21;
    double average, dDisttanceDd12, dDisttanceDd21;

    public PeptideBondDistance(PeptideBond pb1, PeptideBond pb2) {
        this.pb1 = pb1;
        this.pb2 = pb2;
        pb1Hpb2O = new FreeDistance(pb1.H, pb2.O);
        pb2Hpb1O = new FreeDistance(pb2.H, pb1.O);
        hhDistance = new FreeDistance(pb1.H, pb2.H);
        ooDistance = new FreeDistance(pb1.O, pb2.O);
        update();
        //System.out.println("Created "+this);
    }
    /* This method was used to test the derivation of the distance (and was useful in finding bugs
    public void update() {
        double epsilon = 0.0000001;
        update1();
        double saveDistance = distance;
        double saveDerivativeX12 = dDistanceDx12;
        double saveDerivativeY12 = dDistanceDy12;
        double saveDerivativeZ12 = dDistanceDz12;
        double saveDerivativeX21 = dDistanceDx21;
        double saveDerivativeY21 = dDistanceDy21;
        double saveDerivativeZ21 = dDistanceDz21;
        pb1.H.addToX(epsilon);
        update1();
        System.out.println((distance-saveDistance)/epsilon);
        System.out.println("x12 "+saveDerivativeX12);
        System.out.println("y12 "+saveDerivativeY12);
        System.out.println("z12 "+saveDerivativeZ12);
        System.out.println("x21 "+saveDerivativeX21);
         System.out.println("y21 "+saveDerivativeY21);
         System.out.println("z21 "+saveDerivativeZ21);
          if (1 == 1) throw new RuntimeException();                 
    }                 */

    public void update() {
        pb1Hpb2O.update();
        pb2Hpb1O.update();

        d12 = pb1Hpb2O.distance();
        d12ALPHA = d12 * ALPHA;
        expD12ALPHA = Math.exp(d12ALPHA);
        dExpD12ALPHAdD12 = ALPHA * expD12ALPHA;

        d21 = pb2Hpb1O.distance();
        d21ALPHA = d21 * ALPHA;
        expD21ALPHA = Math.exp(d21ALPHA);
        dExpD21ALPHAdD21 = ALPHA * expD21ALPHA;

        average = (expD12ALPHA + expD21ALPHA) / 2;
        distance = Math.log(average) / ALPHA;
        dDisttanceDd12 = 0.5 / ALPHA / average * dExpD12ALPHAdD12;
        dDisttanceDd21 = 0.5 / ALPHA / average * dExpD21ALPHAdD21;
        dDistanceDx12 = dDisttanceDd12 * pb1Hpb2O.dDistanceDx();
        dDistanceDy12 = dDisttanceDd12 * pb1Hpb2O.dDistanceDy();
        dDistanceDz12 = dDisttanceDd12 * pb1Hpb2O.dDistanceDz();
        dDistanceDx21 = dDisttanceDd21 * pb2Hpb1O.dDistanceDx();
        dDistanceDy21 = dDisttanceDd21 * pb2Hpb1O.dDistanceDy();
        dDistanceDz21 = dDisttanceDd21 * pb2Hpb1O.dDistanceDz();

    }

    public double distance() {
        return distance;
    }

    public double dDistanceDx12() {
        return dDistanceDx12;
    }

    public double dDistanceDy12() {
        return dDistanceDy12;
    }

    public double dDistanceDz12() {
        return dDistanceDz12;
    }

    public double dDistanceDx21() {
        return dDistanceDx21;
    }

    public double dDistanceDy21() {
        return dDistanceDy21;
    }

    public double dDistanceDz21() {
        return dDistanceDz21;
    }

    public String toString() {
        return "PeptideBondDistance (" + pb1 + ";" + pb2 + ") " + distance;
    }

}
