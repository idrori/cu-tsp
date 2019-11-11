/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public interface Updateable {
    public void update(int numberOfUpdates) throws UpdateableException;

    public void setNumberOfUpdates(int numberOfUpdates);
}
