/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package programs;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 13/12/2010
 * Time: 14:16:13
 * To change this template use File | Settings | File Templates.
 */
import meshi.dataStructures.Set;
import meshi.dataStructures.SetAsArray;
import meshi.design.*;
public class Design {
    public static void main(String[] argv) {
        ConstantFunction cf = new ConstantFunction(1.000010000);
        LinearFunction lf  = new LinearFunction(2,3) ;
        SinFunction sf = new SinFunction(2,2,Math.PI);
        for(double x = 0 ; x < 2*Math.PI; x = x+0.1) {
            System.out.println(x+"\t"+
                    cf.valueAt(x)+"\t"+cf.derivativeAt(x)+"\t"+
                    lf.valueAt(x)+"\t"+lf.derivativeAt(x)+"\t"+
                    sf.valueAt(x)+"\t"+sf.derivativeAt(x)+"\t");
        }
         System.out.println(Math.pow(2,0.5));
        System.out.println(Math.pow(-2,0.5)) ;
        double zero = 0.0;
        System.out.println(-2/zero) ;
        System.out.println(Double.isNaN(Math.pow(2,0.5)));
       System.out.println(Double.isNaN(Math.pow(-2,0.5))) ;
       System.out.println(Double.isInfinite(-2/zero)) ;

       Set set = new SetAsArray();
        set.add(new Integer(2));
        set.add(new Double(2));
        set.add(new Integer(3));
        set.add(new Integer(3));
        System.out.println(set.size());
        System.out.println(set.contains(new Integer(3)));
        System.out.println(set.contains(new Integer(4)));        

       }
}
