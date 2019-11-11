package meshi.util.filters;

import java.io.File;
import java.io.FileFilter;

/**
 * Created by chen on 03/12/2015.
 */
public class PDBfileFilter  implements FileFilter {
    public boolean accept(File file) {
        return file.getAbsolutePath().endsWith(".pdb");
    }
}
