/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public class AttributesRack {
    public static final int MAX_ATTRIBUTES = 27;
    private final Object[] internalArray = new Object[MAX_ATTRIBUTES];

    public AttributesRack() {
        for (Object o : internalArray) {
            o = null;
        }
    }

    public void addAttribute(MeshiAttribute attribute) {
        int key = attribute.key();
        if (key > MAX_ATTRIBUTES) throw new RuntimeException("Key = " + key + " is larger than " +
                "MAX_ATTRIBUTES " + MAX_ATTRIBUTES);
//	if (internalArray[key] != null) throw new RuntimeException("Cannot add attribute "+attribute+
//								   "Key = "+key+" is aleady in use ");
        internalArray[key] = attribute;
    }

    public MeshiAttribute getAttribute(int key) {
        return (MeshiAttribute) internalArray[key];
    }
}
 
