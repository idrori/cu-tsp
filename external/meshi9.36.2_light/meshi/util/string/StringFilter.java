/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.string;

import meshi.util.*;
import meshi.util.filters.*;

import java.util.*;

public abstract class StringFilter implements Filter {
    protected StringList keys;

    public StringFilter(String key) {
        this(new StringList(key));
    }

    public StringFilter(StringList keys) {
        this.keys = keys;
    }

    public boolean accept(Object obj) {
        if (obj instanceof String) {
            Iterator keysIterator = keys.iterator();
            String key;
            while ((key = (String) keysIterator.next()) != null)
                if (accept((String) obj, key)) return true;
            return false;
        }
        throw new MeshiException("Tried to StringFilter:\n" +
                obj + " of class: " + obj.getClass());
    }

    public abstract boolean accept(String string, String key);
}
