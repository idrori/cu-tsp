package meshi.util.proteinReaders;

import meshi.molecularElements.Protein;
import meshi.molecularElements.ProteinMetaData;

/**
 * Created by chen on 01/12/2015.
 */
public interface ProteinReader {
    public Protein readProtein(ProteinMetaData proteinMetaData);
}
