package meshi.util.filters;

import meshi.molecularElements.atoms.Atom;
import meshi.parameters.ResidueType;

/**
 * Created by chen on 12/12/2015.
 */
public class AromaticSideChainAtom implements Filter {
    public boolean accept(Object obj){
        Atom atom = (Atom) obj;
        if (atom.isBackbone()) return false;
        if (atom.residue().type == ResidueType.PHE) return true;
        if (atom.residue().type == ResidueType.TYR) return true;
        if (atom.residue().type == ResidueType.TRP) return true;
        if (atom.residue().type == ResidueType.HIS) return true;
        return false;
    }
}
