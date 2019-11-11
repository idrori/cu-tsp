/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.string;

import meshi.util.*;
import meshi.util.filters.*;;

public class StringFilterGrep extends StringFilter {
    public StringFilterGrep(String key) {
        super(key);
    }

    public StringFilterGrep(StringList keys) {
        super(keys);
    }

    public boolean accept(String string, String key) {
        return (string.indexOf(key) > -1);
    }
}
	
    
