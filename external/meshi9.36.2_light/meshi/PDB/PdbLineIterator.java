package meshi.PDB;

import meshi.util.file.MeshiLineReader;

import java.io.IOException;
import java.util.Iterator;

/**
 * Created by chen on 15/07/2017.
 */
public class PdbLineIterator implements Iterator <PdbLine> {
    private PdbReader reader;
    private PdbLine nextLine;
    public PdbLineIterator(PdbReader reader){
        this.reader = reader;
        nextLine = reader.readPdbLine();
    }

    public PdbLine next() {
        PdbLine out = nextLine;
            nextLine = reader.readPdbLine();
        return out;
    }
    public boolean hasNext() {
        if (nextLine == null) return false;
        return (!nextLine.isEnd());
    }

    public MeshiLineReader sourceFile() {
        return reader;
    }

    public void remove() {throw new RuntimeException("Unsupported operation");}

}
