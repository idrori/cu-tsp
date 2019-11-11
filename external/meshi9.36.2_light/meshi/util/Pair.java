/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

public class Pair {
    public final String key;
    private Object value;

    public Pair(String key, Object value) {
        this.key = key;
        this.value = value;
    }

    public Pair(String key, double value) {
        this(key, new Double(value));
    }

    public Object value() {
        return value;
    }

    public double doubleValue() {
        if (!(value instanceof Double)) throw new RuntimeException("Cannot convert " + value + "to a double");
        return (Double) value;
    }

    public void addToValue(double number) {
        if (!(value instanceof Double)) throw new RuntimeException("Cannot add number to " + value);
        value = new Double(number + (Double) value);
    }

    public void setValue(Object obj) {
        value = obj;
    }

    public void setValue(double number) {
        value = new Double(number);
    }

    public String toString() {
        if (value instanceof Double)
            return String.format("%-20s%-20.4f", key, value);
        else return String.format("%-20s%-20s", key, value);
    }
}
