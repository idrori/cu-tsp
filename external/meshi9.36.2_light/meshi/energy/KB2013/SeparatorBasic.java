package meshi.energy.KB2013;

import meshi.molecularElements.atoms.Atom;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 05/08/13
 * Time: 08:25
 * To change this template use File | Settings | File Templates.
 */
public class SeparatorBasic extends Separator {
    public int numberOfSeparations() {return 2;}
    public int separation(Atom atom1, Atom atom2) {
        int s = atom1.residue().number()-atom2.residue().number();
        if (s == 0)
            throw new RuntimeException("This is weird.\n"+atom1+"\n"+atom2);
        if (s < 0) s = -s;
        if (s <= 2) return 0;
        return 1;
    }
    public String name() {
        return "SeparatorBasic";
    }
}
