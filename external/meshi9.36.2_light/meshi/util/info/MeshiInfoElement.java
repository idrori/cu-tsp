/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.energy.AbstractEnergy;
import meshi.util.Utils;

/**
 *
 */
public abstract class MeshiInfoElement {
    public final String comment;
    public final InfoType type;
    public static int FIELD_LENGTH = 30;

    public MeshiInfoElement(InfoType type, String comment) {
        this.type = type;
        this.comment = comment;
    }

    public MeshiInfoElement(MeshiInfoElement element) {
        this(element.type, element.comment);
    }

    public InfoValueType valueType() {
        return type.valueType;
    }

    public boolean equals (Object other) {
        return type.equals(((MeshiInfoElement)other).type);
    }

    public String toString() {
        return (type.toString() + "                                                                               ").substring(0, FIELD_LENGTH) + "\t" + comment;
    }

    public double doubleValue() {
        if (type.valueType == InfoValueType.DOUBLE)
            return ((DoubleInfoElement) this).value();
        throw new RuntimeException("No doubleValue for " + type.valueType);
    }

    public int intValue() {
        if (type.valueType == InfoValueType.INTEGER)
            return ((IntInfoElement) this).value();
        throw new RuntimeException("No intValue for " + type.valueType);
    }

    public String stringValue() {
        if (type.valueType == InfoValueType.STRING)
            return ((StringInfoElement) this).value();
        throw new RuntimeException("No stringValue for " + type.valueType);
    }

    public String toXML() {
        String out;
        if (type.name().contains("SsS_PARAM")){
            out = "<" + type.tag.replace('<','o').replace('>','c').replace('+','p').replace('-','m') + "  value=\"";
        }
        else out = "<" + type.tag + "  value=\"";
        //String out = "<" + type.tag + "  value=\"";
        if (type == InfoType.SOLVATION_ENTROPY) {
        }
        switch (type.valueType) {
            case INTEGER:
                out += intValue();
                break;
            case DOUBLE:
                out += doubleValue();
                break;
            case STRING:
                out += stringValue();
                break;
            default:
                throw new RuntimeException("Weird element" + this + " of unknown type " + type.valueType);
        }
        out += "\"/>";
        return out;
    }
}
