/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

/**

 */
public class IntInfoElement extends MeshiInfoElement {
    private int value;

    public IntInfoElement(InfoType type, String comment, int value) {
        super(type, comment);
        value = value;
    }

    public IntInfoElement(IntInfoElement old) {
        this(old.type, old.comment, old.value);
    }


    public void toZero() {
        value = 0;
    }

    public void setValue(int value) {
        this.value = value;
    }

    public int value() {
        return value;
    }

    public String toString() {
        return ((new Integer(value)).toString() + "                                                                               ").substring(0, FIELD_LENGTH) + super.toString();
    }
}
