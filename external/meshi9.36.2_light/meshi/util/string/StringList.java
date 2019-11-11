/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.string;

import meshi.util.*;
import meshi.util.file.*;
import meshi.util.filters.*;

import java.util.*;
import java.io.*;


// constructors
public class StringList extends ArrayList<String> {
    private String comment = null;

    public StringList() {
        super();
    }

    public StringList(String string) {
        super();
        add(string);
    }

    public StringList(StringTokenizer st) {
        super();
        while (st.hasMoreTokens()) add(st.nextToken());
    }

    public StringList(MeshiLineReader MLR) {
        super();
        setComment(MLR.path());
        String string;
        try {
            while ((string = MLR.readLine()) != null)
                add(string);
            MLR.close();
        }
        catch (Exception e) {
            throw new MeshiException("StringList(MeshiLineReader MLR) Error:\n" +
                    e.getMessage());
        }
    }

    public StringList(File file) {
        this(new MeshiLineReader(file));
    }

    public StringList(MeshiLineReader MLR, Filter filter) {
        super();
        setComment(MLR.path());
        String string;
        try {
            while ((string = MLR.readLine()) != null)
                if (filter.accept(string)) add(string);
        }
        catch (Exception e) {
            throw new MeshiException("StringList(MeshiLineReader MLR) Error:\n" +
                    e.getMessage());
        }
    }


    public String lastString() {
        return get(size() - 1);
    }

    public StringList filter(Filter filter) {
        StringList out = new StringList();
        for (String s : this)
            if (filter.accept(s)) out.add(s);
        return out;
    }


    public String startsWith(String key) {
        for (String string : this) {
            if (string.startsWith(key)) return string;
        }
        return null;
    }

    public StringList filterStartsWith(String key) {
        StringFilterStartsWith filter =
                new StringFilterStartsWith(key);
        return filter(filter);
    }

    public StringList filterStartsWith(StringList keys) {
        StringFilterStartsWith filter =
                new StringFilterStartsWith(keys);
        return filter(filter);
    }

    public StringList filterEndsWith(String key) {
        StringFilterEndsWith filter =
                new StringFilterEndsWith(key);
        return filter(filter);
    }

    public StringList filterEndsWith(StringList keys) {
        StringFilterEndsWith filter =
                new StringFilterEndsWith(keys);
        return filter(filter);
    }

    public StringList filterGrep(String key) {
        StringFilterGrep filter =
                new StringFilterGrep(key);
        return filter(filter);
    }

    public StringList filterGrep(StringList keys) {
        StringFilterGrep filter =
                new StringFilterGrep(keys);
        return filter(filter);
    }

    public StringList stringParseAt(int index, StringList separators) {
        return StringParser.bySeparators(get(index), separators);
    }

    public StringList stringParseAt(int index, String separator) {
        return StringParser.bySeparator(get(index), separator);
    }

    public StringList stringParseAt(int index) {
        return StringParser.standard(get(index));
    }
//     public StringList flatten() {
// 	StringList newList = new StringList(),tempList;
// 	Iterator SI = iterator();
// 	Iterator SI1;
// 	String string;
// 	while ((tempList = StringParser.standard((String) SI.next())) != null)
// 	    {
// 		SI1 = tempList.iterator();
// 		while ((string = (String) SI1.next()) != null)
// 		    newList.add(string);
// 	    }
// 	return newList;
//     }


    public static StringList standardSeparators() {
        StringList newList = new StringList();
        newList.add(" ");
        newList.add(",");
        newList.add(";");
        newList.add("\n");
        newList.add("\t");
        return newList;
    }

    static class IsString implements Filter {
        public boolean accept(Object obj) {
            return (obj instanceof String);
        }
    }

    public void setComment(String s) {
        comment = s;
    }

    public void print(MeshiWriter writer) {
        for (String s : this)
            writer.println(s);
    }
}	
		
	
