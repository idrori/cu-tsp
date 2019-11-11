package meshi.molecularElements;

import meshi.molecularElements.extendedAtoms.Pro;
import meshi.util.Utils;
import meshi.util.info.MeshiInfo;

import java.util.*;

/**
 * Created by chen on 01/12/2015.
 */
public class ProteinMetaData extends Hashtable<ProteinMetaData.MetaDataKey,Object> {
    public enum MetaDataKey {NAME, FILE_NAME, DIRECTORY, TARGET}
    String fileName;
    String directory;
    String target;

    public ProteinMetaData() {
        super();
    }

    public ProteinMetaData(String directory, String fileName, String target) {
        put(MetaDataKey.FILE_NAME, fileName);
        put(MetaDataKey.DIRECTORY, directory);
        put(MetaDataKey.TARGET, target);
    }

    public String toString() {
        String out = "ProteinMetaData ";
        for (MetaDataKey key : MetaDataKey.values())
            out += key.toString()+"="+get(key)+" ";
        out += "\n";
        return out;
    }
    public ProteinMetaData duplicate() {
        ProteinMetaData out = new ProteinMetaData();
        for (MetaDataKey key : MetaDataKey.values()) {
            Object value = get(key);
            if (value != null)
                out.put(key,value);
        }
        return out;
    }

}
