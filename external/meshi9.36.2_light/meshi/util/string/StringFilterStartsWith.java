/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.string;

import meshi.util.*;
import meshi.util.filters.*;;

public class StringFilterStartsWith extends StringFilter {
    public StringFilterStartsWith(String key) {
        super(key);
    }

    public StringFilterStartsWith(StringList keys) {
        super(keys);
    }

    public boolean accept(String string, String key) {
        return string.startsWith(key);
    }
}
	
    
