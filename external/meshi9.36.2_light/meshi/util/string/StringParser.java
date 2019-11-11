/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.string;

import meshi.util.*;

import java.util.*;

public class StringParser {
    public static StringList bySeparators(String string, StringList separators) {
        StringList newList = new StringList();
        int ptr, minPtr;
        String minSep = "";
        String parseMe;
        if (string == null) return null;
        parseMe = string;
        while (parseMe.length() > 0) {
            minPtr = 10000;
            for (String separator : separators) {
                ptr = parseMe.indexOf(separator);
                if (ptr != -1)
                    if (ptr < minPtr) {
                        minPtr = ptr;
                        minSep = separator;
                    }
            }
            if (minPtr == 10000) {
                newList.add(parseMe);
                parseMe = "";
            } else if (minPtr == 0) {
                parseMe = parseMe.substring(minSep.length());
            } else {
                newList.add(parseMe.substring(0, minPtr));
                parseMe = parseMe.substring(minPtr + minSep.length());
            }
        }
        return newList;
    }

    public static StringList bySeparator(String string, String separator) {
        return bySeparators(string, new StringList(separator));
    }

    public static StringList standard(String string) {
        return bySeparators(string, StringList.standardSeparators());
    }

    public static StringList byPosition(String string, ArrayList<Integer> positions) {
        StringList newList = new StringList();
        int from = 0, to;
        String parseMe;
        parseMe = string;
        for (Integer position : positions) {
            to = position.intValue();
            if (to + 1 > parseMe.length()) {
                newList.add(parseMe.substring(from, parseMe.length()));
                return newList;
            } else {
                newList.add(parseMe.substring(from, to + 1));
                from = to + 1;
            }
        }
        newList.add(parseMe.substring(from, parseMe.length()));
        return newList;
    }

    public static StringList breakPath(String path) {
        return bySeparator(path, "/");
    }

    public static StringList breakFileName(String fileName) {
        return bySeparator(fileName, ".");
    }
}
