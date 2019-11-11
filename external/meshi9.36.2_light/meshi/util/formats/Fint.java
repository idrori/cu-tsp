/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.formats;

import meshi.util.*;

import java.util.*;

public class Fint extends Format {
    public static final Fint STANDARD = new Fint(12, -1);
    public static final Fint SHORT = new Fint(8, -1);
    public static final Fint SHORTER = new Fint(7, -1);

    public Fint(int field, int indent) {
        super(field, indent);
    }

    public String f(Object value) {
        return f(((Integer) value).intValue());
    }

    public String f(int value) {
        if (INDENT == -1) return fintL(value, FIELD);
        else if (INDENT == 1) return fintR(value, FIELD);
        else throw new MeshiException("unknown indentation flag " + INDENT);
    }
}
    
