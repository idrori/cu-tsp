package programs;

import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.*;
import java.net.URL;

/**
 * Created by chen on 29/11/2015.
 */
public class GetStructuresFromPDB {
    public static void main(String[] args) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(args[0]);
        String line;
        reader.readLine();
        while ((line = reader.readLine()) != null) {
            String[] words = line.split(" ");
            String proteinName = words[0].substring(0, 4) + ".pdb";
            File file = new File(proteinName);
            if (!file.exists()) {
                URL url = new URL(" \thttp://www.rcsb.org/pdb/files/" + proteinName);
                System.out.println(url.toString());
                InputStream inputStream = url.openStream();
                InputStreamReader inputStreamReader = new InputStreamReader(inputStream);
                BufferedReader bufferedReader = new BufferedReader(inputStreamReader);
                MeshiWriter writer = new MeshiWriter(proteinName);
                while ((line = bufferedReader.readLine()) != null)
                    writer.println(line);
                writer.close();
                bufferedReader.close();
                inputStreamReader.close();
                inputStream.close();
            }
        }
    }
}
