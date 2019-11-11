/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.formats;

import meshi.util.*;

import java.util.*;

public class Fdouble extends Format {
    public final int PRECISION;
    public static final Fdouble STANDARD = new Fdouble(8, 12, -1);
    public static final Fdouble SHORT = new Fdouble(4, 8, -1);
    public static final Fdouble SHORTER = new Fdouble(4, 7, -1);
    public static final Fdouble VERY_SHORT = new Fdouble(4, 4, -1);

    public Fdouble(int precision, int field, int indent) {
        super(field, indent);
        PRECISION = precision;
    }

    public String f(Object value) {
        return f(((Double) value).doubleValue());
    }

    public String f(double value) {
        if (INDENT == -1) return fdoubleL(value, PRECISION, FIELD);
        else if (INDENT == 1) return fdoubleR(value, PRECISION, FIELD);
        else throw new MeshiException("unknown indentation flag " + INDENT);
    }
}
    
