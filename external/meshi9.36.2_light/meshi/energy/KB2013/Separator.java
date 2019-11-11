package meshi.energy.KB2013;

import meshi.molecularElements.atoms.Atom;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 05/08/13
 * Time: 00:46
 * To change this template use File | Settings | File Templates.
 */
public abstract class Separator {
    public abstract int numberOfSeparations();
    public abstract int separation(Atom atom1,Atom atom2);
    public abstract String name();
    public static Separator getSeparator(String name) {
        if (name.equals((new SeparatorBasic()).name()))
            return new SeparatorBasic();
        if (name.equals((new Separator10()).name()))
                    return new Separator10();
        throw new RuntimeException(name+" is not a name of a valid Separator");
    }
}
