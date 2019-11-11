/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

/**
 * Created by IntelliJ IDEA.
 * User: keasar
 * Date: 06/11/2010
 * Time: 10:58:29
 * To change this template use File | Settings | File Templates.
 */
public enum InfoValueType {
    DOUBLE, INTEGER, STRING, OTHER;


    public static boolean numeric(InfoValueType type) {
        return ((type == DOUBLE) || (type == INTEGER));
    }

    public static boolean numeric(MeshiInfoElement element) {
        return numeric(element.type.valueType);
    }

}
