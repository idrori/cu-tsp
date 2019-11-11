/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.util.Utils;

import java.util.ArrayList;

/**
 *
 */
public class MeshiInfo {
    public static final int FIELD_LENGTH = 30;

    public final String comment;
    public final InfoType type;
    private Object value;

    public ArrayList<MeshiInfo> getChildren() {
        return children;
    }


    private final ArrayList<MeshiInfo> children;

    public MeshiInfo(InfoType type, Object value, String comment, MeshiInfo[] children) {
        this.children = new ArrayList<MeshiInfo>();
        this.type = type;
        this.value = value;
        this.comment = comment;
        if (children != null) {
            for (MeshiInfo child : children) {
                this.children.add(child);
            }
        }
    }

    public MeshiInfo(InfoType type, Object value,  String comment) {
        this(type, value, comment, null);
    }
    public MeshiInfo(InfoType type, String comment) {
        this(type, null, comment, null);
    }

    public ArrayList<MeshiInfo> flatten() {
        ArrayList<MeshiInfo> flattenned = new ArrayList();
        flatten(this, flattenned);
        return flattenned;
    }

    private static void flatten(MeshiInfo info, ArrayList<MeshiInfo> flattenned){
        flattenned.add(new MeshiInfo(info.type, info.value, info.comment));
        if (info.children != null)
            for (MeshiInfo child : info.children)
                if (child != null)
                    flatten(child,flattenned);
    }

    public Object getValue() {
        return value;
    }

    public void setValue(Object value) {
        this.value = value;
    }

    public void reset() {
        setValue(new Double(0.0));
        for (MeshiInfo child : children)
            child.reset();
    }
    public String toString() {
        return "MeshiInfo " + type + " " + getValue();
    }



    public String toXml(String indentation) {
        String out;
        if (type != null)
            out = indentation + "<" + type + " value=\"" + getValue() + "\"/>";
        else  out = indentation + "<" + comment + " value=\"" + getValue() + "\"/>";
        if (children != null)
            for (MeshiInfo info : children)
                out += info.toXml(indentation);
        return out;
    }
    public String toXml() {
        return toXml("");
    }
}
