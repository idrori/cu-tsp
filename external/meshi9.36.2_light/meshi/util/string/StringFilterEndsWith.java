/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.string;

import meshi.util.*;
import meshi.util.filters.*;;

public class StringFilterEndsWith extends StringFilter {
    public StringFilterEndsWith(String key) {
        super(key);
    }

    public StringFilterEndsWith(StringList keys) {
        super(keys);
    }

    public boolean accept(String string, String key) {
        return string.endsWith(key);
    }
}
	
    
