package meshi.molecularElements;

import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Created by chen on 12/09/2016.
 */
public class MultiModelEnsemble extends ArrayList<Protein>{
    private enum Status {IN_MODEL, OUT}
    public MultiModelEnsemble() {}

    public MultiModelEnsemble(String fileName) throws IOException{
        super();
        if (isMultiModelEnsemble(fileName)) {
            MeshiLineReader reader = new MeshiLineReader(fileName);
            String line;
            Status status = Status.OUT;
            MeshiWriter tmpWriter = null;
            String modelName = null;
            int modelNumber = 1;
            int nLine = 10000;
            while ((line = reader.readLine()) != null) {
                if (status == Status.OUT) {
                    if (line.startsWith("MODEL")) {
                        String[] words = line.split("\\s+");
                        modelName = words[1];
                        nLine = 0;
                        System.out.println(line);
                        status = Status.IN_MODEL;
                        tmpWriter = new MeshiWriter("temp.pdb");
                    }
                } else {
                    if (line.equals("ENDMDL")) {
                        if (nLine == 0)
                            throw new RuntimeException("This is weird. Empty Protein. ");
                        status = Status.OUT;
                        if (tmpWriter != null) tmpWriter.close();
                        Protein tmpProtein = new Protein(new AtomList("temp.pdb"), new ResidueExtendedAtomsCreator());
                        if (modelName == null)
                            throw new RuntimeException("this is weird.");
                        tmpProtein.setName(modelName);
                        modelNumber++;
                        add(tmpProtein);
                    } else {
                        tmpWriter.println(line);
                        nLine++;
                    }
                }
            }
            reader.close();
        }
        else add(new Protein(new AtomList(fileName), new ResidueExtendedAtomsCreator()));
    }
    public void print(MeshiWriter writer) {
        for (Protein protein : this) {
          if (protein != null) {
              writer.println("MODEL " + protein.name());
              protein.atoms().print(writer);
              writer.println("ENDMDL");
          }
        }
        writer.println("END");
    }

    public static boolean isMultiModelEnsemble(String fileName) throws IOException{
        boolean out = false;
        MeshiLineReader reader = new MeshiLineReader(fileName);
        String line;
        Status status = Status.OUT;
        while ((line = reader.readLine()) != null) {
            if (status == Status.OUT) {
                if (line.startsWith("MODEL")) {
                    status = Status.IN_MODEL;
                }
            }
            else {
                if (line.startsWith("END")){
                    status = Status.OUT;
                    out = true;
                }
            }
        }
        reader.close();
        return out;
    }
}
