package programs;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

/**
 * Created by IntelliJ IDEA.
 * User: chen
 * Date: 28/08/13
 * Time: 22:46
 * To change this template use File | Settings | File Templates.
 */
public class FindDiscontinuousDomains {
    public static void main(String[] args) throws IOException{
        File here = new File(".");
        File discontinuous = new File("discontinuous");
        File[] files = here.listFiles();
        for (File dir : files) {
            if (dir.isDirectory()){
                File[] proteinFiles = dir.listFiles();
                for (File proteinFile : proteinFiles) {
                    if (proteinFile.getName().endsWith(".N.pdb")){
                        System.out.println(proteinFile.getAbsolutePath());
                        if (GetDomains.discontinuous(proteinFile)) {
                            Files.move(proteinFile.toPath(),discontinuous.toPath());
                            System.out.println(proteinFile.getAbsolutePath()+" discontinuous ");
                        }
                    }
                }
            }
        }

    }

}
