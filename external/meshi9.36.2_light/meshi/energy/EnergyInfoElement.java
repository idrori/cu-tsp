/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy;

import meshi.util.Utils;
import meshi.util.info.*;

/**

 */
public class EnergyInfoElement extends MeshiInfo {
    private double weight = Double.MIN_VALUE;
    public double weight() {
        return weight;
    }


    public EnergyInfoElement(InfoType type, String comment, double weight, MeshiInfo[] infoList) {
        this(type, comment, weight,infoList,0);
    }
    public EnergyInfoElement(InfoType type, String comment, double weight) {
        super(type, null, comment, null);
        this.weight = weight;
    }
    public EnergyInfoElement(InfoType type) {
        this(type, type.tag);
        setValue(new Double(0));
    }
    public EnergyInfoElement(InfoType type, String comment) {
        super(type, null, comment, null);
        this.weight = 1;
    }
    public EnergyInfoElement(EnergyInfoElement energyInfoElement) {
        this(energyInfoElement.type,
                energyInfoElement.comment,
                energyInfoElement.weight,
                (MeshiInfo[]) energyInfoElement.getChildren().toArray(),
                ((Double)energyInfoElement.getValue()).doubleValue());
    }
    public EnergyInfoElement(InfoType type,
                             String comment,
                             double weight,
                             MeshiInfo[] children,
                             double value) {
        super(type, Double.valueOf(value), comment, children);
        this.weight = weight;
    }


    public double energy() {
        return ((Double) getValue()).doubleValue();
    }

}
