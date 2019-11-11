package programs;

import meshi.energy.EnergyInfoElement;
import meshi.energy.compatebility.StaticFeaturesCreator;
import meshi.energy.compatebility.StaticFeatures;
import meshi.molecularElements.Protein;
import meshi.molecularElements.extendedAtoms.ResidueExtendedAtomsCreator;
import meshi.util.CommandList;
import meshi.util.MeshiProgram;
import meshi.util.Utils;
import meshi.util.info.MeshiInfo;

import java.io.File;
import java.io.IOException;

/**
 * Created by chen on 23/03/2016.
 */
public class StaticFeaturesTest extends MeshiProgram{
    public static void main(String[] argv) throws IOException{
        initRandom(0);
        CommandList commands = new CommandList(argv[0]);
        File file = new File(argv[1]);
        Protein model = Utils.getProtein(commands, file.getAbsolutePath(), new ResidueExtendedAtomsCreator(), Utils.defaultExceptionHandler);
        Utils.AssignDSSP(model,argv[2]);
        StaticFeaturesCreator creator = new StaticFeaturesCreator();
        StaticFeatures staticFeatures = (StaticFeatures) creator.createEnergyTerm(model,null,commands) ;
        EnergyInfoElement info = staticFeatures.evaluate();
        System.out.println(model.name()+" "+info.getValue());
        for (MeshiInfo meshiInfo : info.getChildren()) {
            System.out.println(meshiInfo.comment+" "+meshiInfo.getValue());
        }
    }
}
