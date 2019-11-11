package meshi.energy.KB2013;

import meshi.molecularElements.atoms.Atom;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 05/08/13
 * Time: 08:25
 * To change this template use File | Settings | File Templates.
 */
public class Separator10 extends Separator {
    public int numberOfSeparations() {return 20;}
    public int separation(Atom atom1, Atom atom2) {
        int s = atom1.residue().number()-atom2.residue().number();
        if (s == 0) throw new RuntimeException("This is weird.");
        if (s < 0) s = -s;
        if (s <= 10) return s-1;
        if (s >= 101)return 19;
        return (int) Math.ceil(s/10.0)+8;
    }
    public String name() {
        return "Separator10";
    }


}
