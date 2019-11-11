package meshi.molecularElements.loops;

import meshi.molecularElements.Residue;

/**
* Created by chen on 12/02/2015.
*/
public class Three {
    public Residue getFirst() {
        return first;
    }

    private Residue first = null;

    public Residue getSecond() {
        return second;
    }

    private Residue second = null;

    public Residue getThird() {
        return third;
    }

    public void setThird(Residue third) {
        this.third = third;
    }

    public void setSecond(Residue second) {
        this.second = second;
    }

    public void setFirst(Residue first) {
        this.first = first;
    }

    private Residue third = null;
public String toString() {return "Three residues: "+first+" "+second+" "+third;}
}
