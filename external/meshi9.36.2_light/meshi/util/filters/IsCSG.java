package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.AtomType;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 14/06/12
 * Time: 15:27
 * To change this template use File | Settings | File Templates.
 */
public class IsCSG implements Filter{
    public static final IsCSG isCSG = new IsCSG();
    public boolean accept(Object obj){
        if (!(obj instanceof Atom))return false;
        return ((Atom)obj).type() == AtomType.CSG;
    }
}
