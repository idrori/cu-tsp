package meshi.util.dssp;

import meshi.molecularElements.Protein;
import meshi.util.dssp.DSSP;

/**
 * Created by chen on 14/12/2015.
 */
public interface DsspReader {
    public DSSP read(Protein protein);
}
