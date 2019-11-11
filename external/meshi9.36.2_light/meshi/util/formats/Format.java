/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.formats;

import meshi.util.*;

public class Format {
    public final int FIELD;
    public final int INDENT;
    public static final Format STANDARD = new Format(12, -1);
    public static final Format SHORTER = new Format(7, -1);
    public static final Format SHORT = new Format(8, -1);
    public static final Format VERY_SHORT = new Format(4, -1);

    public Format(int field, int indent) {
        FIELD = field;
        INDENT = indent;
    }

    public String f(Object value) {
        if (INDENT == -1) return fstringL(value.toString(), FIELD);
        else if (INDENT == 1) return fstringR(value.toString(), FIELD);
        else throw new MeshiException("unknown indentation flag " + INDENT);
    }

    private static String fdouble(double d, int precision, int fieldSize) {
        int i;
        String stringd = (new Double(d)).toString();
        if (!(stringd.charAt(0) == '-')) stringd = " " + stringd;
        int dotIndex = -1;
        int expIndex = -1;
        for (i = 0; i < stringd.length(); i++)
            if (stringd.charAt(i) == 'E') expIndex = i;
        for (i = 0; i < stringd.length(); i++)
            if (stringd.charAt(i) == '.') dotIndex = i;
        if (dotIndex == -1) {
            stringd += ".";
            for (i = 0; i < precision; i++)
                stringd += "0";
            return stringd;
        }
        if (dotIndex > fieldSize) {
            stringd = "";
            for (i = 0; i < fieldSize; i++)
                stringd += "*";
            return stringd;
        }
        if (expIndex == -1) {
            int diff = stringd.length() - dotIndex;
            if (diff > precision) return stringd.substring(0, dotIndex + precision + 1);
            for (i = 0; i < precision - diff + 1; i++)
                stringd += "0";
            return stringd;
        }
        int expDiff = stringd.length() - expIndex;
        int dotDiff = expIndex - dotIndex;
        if (dotDiff < precision) return stringd;
        stringd = stringd.substring(0, dotIndex + precision + 2 - expDiff) +
                stringd.substring(expIndex, stringd.length());
        return stringd;
    }

    public static String fdoubleL(double d, int precision, int fieldSize) {
        return fstringL(fdouble(d, precision, fieldSize), fieldSize);
    }

    public static String fdoubleR(double d, int precision, int fieldSize) {
        return fstringR(fdouble(d, precision, fieldSize), fieldSize);
    }

    private static String fint(int ii, int fieldSize) {
        int i, n;
        String temp = (new Integer(ii)).toString();
        if (temp.length() > fieldSize) {
            temp = "";
            for (i = 0; i < fieldSize; i++)
                temp += "*";
        }
        return temp;
    }

    public static String fintL(int ii, int fieldSize) {
        return fstringL(fint(ii, fieldSize), fieldSize);
    }

    public static String fintR(int ii, int fieldSize) {
        return fstringR(fint(ii, fieldSize), fieldSize);
    }

    private static String fstring(String s, int fieldSize) {
        if (s.length() < fieldSize) return s;
        return s.substring(0, fieldSize);
    }

    public static String fstringL(String s, int fieldSize) {
        String temp = fstring(s, fieldSize);
        for (int i = temp.length(); i < fieldSize; i++)
            temp += " ";
        return temp;
    }

    public static String fstringR(String s, int fieldSize) {
        String temp = fstring(s, fieldSize);
        for (int i = temp.length(); i < fieldSize; i++)
            temp = " " + temp;
        return temp;
    }
}
