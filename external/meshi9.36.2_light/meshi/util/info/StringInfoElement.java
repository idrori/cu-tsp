/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

/**

 */
public class StringInfoElement extends MeshiInfoElement {
    private String value;

    public StringInfoElement(InfoType infoType, String comment, String value) {
        super(infoType, comment);
        this.value = value;
    }

    public StringInfoElement(StringInfoElement element) {
        this(element.type, element.comment, element.value);
    }

    public String value() {
        return value;
    }
}
