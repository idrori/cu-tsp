package meshi.util.info;

import meshi.molecularElements.Residue;
import meshi.util.Utils;


import java.util.ArrayList;

/**
 * Created by chen on 17/08/2017.
 */
public class ResidueInfo extends MeshiInfoElement {
    private Residue residue;
    private ArrayList<DoubleInfoElement> data;

    public ResidueInfo(Residue residue) {
        super(InfoType.ResidueInfo,"Residue data");
        this.residue = residue;
        data = new ArrayList<DoubleInfoElement>();
    }

    public void add(DoubleInfoElement element) {
        data.add(element);
    }
    public String toString() {
        String out = "<ResidueInfo ";
        out += "chain=\""+residue.chain()+"\" ";
        out += "number="+ Utils.F(residue.number());
        out += "type="+Utils.F(residue.type.toString());
        out += "typeID="+Utils.F(residue.type.ordinal());
        if (!residue.dummy()) {
            out += "caX=" + Utils.F(residue.ca().x());
            out += "caY=" + Utils.F(residue.ca().y());
            out += "caZ=" + Utils.F(residue.ca().z());
        }
        for (DoubleInfoElement element : data) {
            boolean toPrint = true;
            if (element.type.equals(InfoType.RESIDUE_CONTACTS8_MCC) & (element.value() == -1))
                toPrint = false;
            if (toPrint)
                out += element.type.tag + "=" + Utils.F(element.value());
        }
        out += "/>";
        return out;
    }

    public boolean dummy() {
        return residue.dummy();
    }
    public int chainNumber() {return residue.getChainNumber();}
    public int number() {return residue.number();}
}
