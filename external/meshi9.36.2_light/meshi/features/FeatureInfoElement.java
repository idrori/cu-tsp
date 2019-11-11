/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.features;

import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfoElement;

import java.util.ArrayList;

/**

 */
public class FeatureInfoElement extends DoubleInfoElement {
    private ArrayList<DoubleInfoElement> infoList;
    private double weight = Double.MIN_VALUE;
    public double weight() {
        return weight;
    }

    public FeatureInfoElement(InfoType type, String comment) {
        this(type, comment, 0,new ArrayList<DoubleInfoElement>());
    }

    public FeatureInfoElement(InfoType type, String comment, double weight) {
        this(type, comment, weight,new ArrayList<DoubleInfoElement>());
    }

    public FeatureInfoElement(InfoType type, String comment, double weight, ArrayList<DoubleInfoElement> infoList) {
        this(type, comment, weight,infoList,0);
    }
    public FeatureInfoElement(InfoType type,
                              String comment,
                              double weight,
                              ArrayList<DoubleInfoElement> infoList,
                              double value) {
        super(type, comment, value);
        this.weight = weight;
        this.infoList = infoList;
    }

    public FeatureInfoElement(FeatureInfoElement other) {
        this(other.type, other.comment, other.weight(),new ArrayList<DoubleInfoElement>(),other.value());
        if (other.infoList()!= null)  {
              for(DoubleInfoElement element : other.infoList())
                  infoList().add(element);
        }
    }

    public double energy() {
        return value();
    }

    public ArrayList<DoubleInfoElement> infoList() {
        return infoList;
    }

    public void resetEnergy() {
        setValue(0);
        if (infoList() != null) {
            for (MeshiInfoElement element : infoList()) {
                ((DoubleInfoElement) element).setValue(0);
            }
        }
    }
}
