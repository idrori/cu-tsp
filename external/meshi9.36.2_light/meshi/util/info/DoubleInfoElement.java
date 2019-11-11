/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyInfoElement;

/**
 *
 */
public class DoubleInfoElement extends MeshiInfoElement {
    private double value;

    public DoubleInfoElement(InfoType type, String comment, double value) {
        super(type, comment);
        this.value = value;
    }

    public DoubleInfoElement(InfoType type, String comment) {
        super(type, comment);
        this.value = Double.MIN_VALUE;
    }

    public DoubleInfoElement(String typeString, String comment, double value) {
        this(InfoType.getType(typeString), comment, value);
    }

    public DoubleInfoElement(String typeString, String comment) {
        this(InfoType.getType(typeString), comment, Double.MIN_VALUE);
    }

    public DoubleInfoElement(DoubleInfoElement old) {
        this(old.type, old.comment, old.value);
    }

    public DoubleInfoElement(EnergyInfoElement old) {
        this(old.type,old.comment, ((Double) old.getValue()).doubleValue());
    }

    public DoubleInfoElement(MeshiInfo old) {
        this(old.type,old.comment, ((Double) old.getValue()).doubleValue());
    }

    public void divideBy(double div) {
        value = value / div;
    }

    public boolean equals(Object other) {
        if (super.equals(other))
            return value == ((DoubleInfoElement) other).value;
        return false;
    }
    public void square() {
        value *= value;
    }

    public void sqrt() {
        value = Math.sqrt(value);
    }

    public void addValue(DoubleInfoElement element) {
        value += element.value;
    }

    public void addValue2(DoubleInfoElement element) {
        value += element.value * element.value;
    }

    public void subtract(DoubleInfoElement element) {
        value -= element.value;
    }

    public void divideBy(DoubleInfoElement element) {
        value /= element.value;
    }

    public void divideBy(IntInfoElement element) {
        value /= element.value();
    }

    public void toZero() {
        value = 0;
    }

    public double setValue(double value) {
        this.value = value;
        return value;
    }

    public double value() {
        return value;
    }

    public String toString() {
        return ((new Double(value)).toString() + "                                                                               ").substring(0, FIELD_LENGTH) + super.toString();
    }


    public static DoubleInfoElement[] generateArray(AbstractEnergy energy) {
        return generateArray(energy, null);
    }

    public static DoubleInfoElement[] generateArray(AbstractEnergy energy, String[] elementNames) {
        int length = 1;
        if (elementNames != null) length += elementNames.length;
        DoubleInfoElement[] out = new DoubleInfoElement[length];
        out[0] = new DoubleInfoElement(InfoType.getType(energy.comment()), energy.comment() + " energy ", -999999);
        for (int i = 1; i < out.length; i++)
            out[i] = new DoubleInfoElement(InfoType.getType(elementNames[i - 1]), energy.comment() + " energy ", -999999);
        return out;
    }
}

